#include "robotoc/parnmpc/unconstr_backward_correction.hpp"

#include <omp.h>
#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

UnconstrBackwardCorrection::UnconstrBackwardCorrection(const OCP& ocp, 
                                                       const int nthreads)
  : nthreads_(nthreads),
    ocp_(ocp),
    intermediate_stage_(ocp.robot, ocp.cost, ocp.constraints),
    terminal_stage_(ocp.robot, ocp.cost, ocp.constraints),
    data_(ocp.N, UnconstrOCPData()),
    performance_index_(),
    corrector_(ocp.N, UnconstrSplitBackwardCorrection(ocp.robot)),
    s_new_(ocp.N+1, SplitSolution(ocp.robot)),
    aux_mat_(ocp.N, Eigen::MatrixXd::Zero(2*ocp.robot.dimv(), 2*ocp.robot.dimv())),
    primal_step_sizes_(Eigen::VectorXd::Zero(ocp.N)),
    dual_step_sizes_(Eigen::VectorXd::Zero(ocp.N)) {
  if (nthreads <= 0) {
    throw std::out_of_range("[UnconstrBackwardCorrection] invalid argument: 'nthreads' must be positive!");
  }
  for (int i=0; i<ocp.N-1; ++i) {
    data_[i] = intermediate_stage_.createData(ocp.robot);
  }
  data_[ocp.N-1] = terminal_stage_.createData(ocp.robot);
}


UnconstrBackwardCorrection::UnconstrBackwardCorrection()
  : nthreads_(0),
    ocp_(),
    intermediate_stage_(),
    terminal_stage_(),
    data_(),
    performance_index_(),
    corrector_(),
    s_new_(),
    aux_mat_(),
    primal_step_sizes_(),
    dual_step_sizes_() {
}


void UnconstrBackwardCorrection::initAuxMat(
    aligned_vector<Robot>& robots, 
    const std::vector<GridInfo>& time_discretization, const Solution& s, 
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  const int N = time_discretization.size() - 1;
  terminal_stage_.evalTerminalCostHessian(robots[0], time_discretization[N], 
                                          s[N-1], data_[N-1], kkt_matrix[N-1], 
                                          kkt_residual[N-1]);
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N; ++i) {
    aux_mat_[i] = kkt_matrix[N-1].Qxx;
  }
  kkt_matrix[N-1].setZero();
  kkt_residual[N-1].setZero();
}


void UnconstrBackwardCorrection::initConstraints(
    aligned_vector<Robot>& robots, 
    const std::vector<GridInfo>& time_discretization, const Solution& s) {
  const int N = time_discretization.size() - 1;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N; ++i) {
    if (i < N-1) {
      intermediate_stage_.initConstraints(robots[omp_get_thread_num()], 
                                          time_discretization[i+1], s[i], data_[i]);
    }
    else {
      terminal_stage_.initConstraints(robots[omp_get_thread_num()], 
                                      time_discretization[i+1], s[i], data_[i]);
    }
  }
}


void UnconstrBackwardCorrection::evalOCP(
    aligned_vector<Robot>& robots, 
    const std::vector<GridInfo>& time_discretization, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
    KKTResidual& kkt_residual) {
  const int N = time_discretization.size() - 1;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N; ++i) {
    if (i == 0) {
      intermediate_stage_.evalOCP(robots[omp_get_thread_num()], 
                                  time_discretization[i+1], q, v, s[i],  
                                  data_[i], kkt_residual[i]);
    }
    else if (i < N-1) {
      intermediate_stage_.evalOCP(robots[omp_get_thread_num()], 
                                  time_discretization[i+1], s[i-1].q, s[i-1].v, 
                                  s[i], data_[i], kkt_residual[i]);
    }
    else {
      terminal_stage_.evalOCP(robots[omp_get_thread_num()], 
                              time_discretization[i+1], s[i-1].q, s[i-1].v, 
                              s[i], data_[i], kkt_residual[i]);
    }
  }
  performance_index_.setZero();
  for (int i=0; i<N; ++i) {
    performance_index_ += data_[i].performance_index;
  }
}


void UnconstrBackwardCorrection::evalKKT(
    aligned_vector<Robot>& robots, 
    const std::vector<GridInfo>& time_discretization, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  const int N = time_discretization.size() - 1;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N; ++i) {
    if (i == 0) {
      intermediate_stage_.evalKKT(robots[omp_get_thread_num()], 
                                  time_discretization[i+1], q, v, s[i], s[i+1], 
                                  data_[i], kkt_matrix[i], kkt_residual[i]);
    }
    else if (i < N-1) {
      intermediate_stage_.evalKKT(robots[omp_get_thread_num()], 
                                  time_discretization[i+1], s[i-1].q, s[i-1].v, 
                                  s[i], s[i+1], data_[i], kkt_matrix[i], kkt_residual[i]);
    }
    else {
      terminal_stage_.evalKKT(robots[omp_get_thread_num()], 
                              time_discretization[i+1], s[i-1].q, s[i-1].v, 
                              s[i], data_[i], kkt_matrix[i], kkt_residual[i]);
    }
  }
  performance_index_.setZero();
  for (int i=0; i<N; ++i) {
    performance_index_ += data_[i].performance_index;
  }
}


void UnconstrBackwardCorrection::coarseUpdate(
    aligned_vector<Robot>& robots, 
    const std::vector<GridInfo>& time_discretization, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  const int N = time_discretization.size() - 1;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N; ++i) {
    if (i == 0) {
      intermediate_stage_.evalKKT(robots[omp_get_thread_num()], 
                                  time_discretization[i+1], q, v, s[i], s[i+1], 
                                  data_[i], kkt_matrix[i], kkt_residual[i]);
      corrector_[i].coarseUpdate(aux_mat_[i+1], time_discretization[i+1].dt, 
                                 kkt_matrix[i], kkt_residual[i], s[i], s_new_[i]);
    }
    else if (i < N-1) {
      intermediate_stage_.evalKKT(robots[omp_get_thread_num()], 
                                  time_discretization[i+1], s[i-1].q, s[i-1].v, 
                                  s[i], s[i+1], data_[i], kkt_matrix[i], kkt_residual[i]);
      corrector_[i].coarseUpdate(aux_mat_[i+1], time_discretization[i+1].dt, 
                                 kkt_matrix[i], kkt_residual[i], s[i], s_new_[i]);
    }
    else {
      terminal_stage_.evalKKT(robots[omp_get_thread_num()], 
                              time_discretization[i+1], s[i-1].q, s[i-1].v, 
                              s[i], data_[i], kkt_matrix[i], kkt_residual[i]);
      corrector_[i].coarseUpdate(time_discretization[i+1].dt, kkt_matrix[i], 
                                 kkt_residual[i], s[i], s_new_[i]);
    }
  }
  performance_index_.setZero();
  for (int i=0; i<N; ++i) {
    performance_index_ += data_[i].performance_index;
  }
}


void UnconstrBackwardCorrection::backwardCorrection(
    const std::vector<GridInfo>& time_discretization, const Solution& s, 
    const KKTMatrix& kkt_matrix, const KKTResidual& kkt_residual, 
    Direction& d) {
  const int N = time_discretization.size() - 1;
  for (int i=N-2; i>=0; --i) {
    corrector_[i].backwardCorrectionSerial(s[i+1], s_new_[i+1], s_new_[i]);
  }
  #pragma omp parallel for num_threads(nthreads_) 
  for (int i=N-2; i>=0; --i) {
    corrector_[i].backwardCorrectionParallel(s_new_[i]);
  }
  for (int i=1; i<N; ++i) {
    corrector_[i].forwardCorrectionSerial(s[i-1], s_new_[i-1], s_new_[i]);
  }
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N; ++i) {
    if (i > 0) {
      corrector_[i].forwardCorrectionParallel(s_new_[i]);
      aux_mat_[i] = - corrector_[i].auxMat();
    }
    UnconstrSplitBackwardCorrection::computeDirection(s[i], s_new_[i], d[i]);
    if (i < N-1) {
      intermediate_stage_.expandPrimalAndDual(time_discretization[i+1].dt,  
                                              kkt_matrix[i], kkt_residual[i], 
                                              data_[i], d[i]);
      primal_step_sizes_.coeffRef(i) = intermediate_stage_.maxPrimalStepSize(data_[i]);
      dual_step_sizes_.coeffRef(i)  = intermediate_stage_.maxDualStepSize(data_[i]);
    }
    else {
      terminal_stage_.expandPrimalAndDual(time_discretization[i+1].dt,  
                                          kkt_matrix[i], kkt_residual[i], 
                                          data_[i], d[i]);
      primal_step_sizes_.coeffRef(i) = terminal_stage_.maxPrimalStepSize(data_[i]);
      dual_step_sizes_.coeffRef(i)  = terminal_stage_.maxDualStepSize(data_[i]);
    }
  }
}


double UnconstrBackwardCorrection::primalStepSize() const {
  return primal_step_sizes_.minCoeff();
}


double UnconstrBackwardCorrection::dualStepSize() const {
  return dual_step_sizes_.minCoeff();
}

void UnconstrBackwardCorrection::integrateSolution(
    const aligned_vector<Robot>& robots,
    const std::vector<GridInfo>& time_discretization, 
    const double primal_step_size, const double dual_step_size,
    Direction& d, Solution& s) {
  assert(robots.size() >= nthreads_);
  const int N = time_discretization.size() - 1;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N; ++i) {
    if (i < N-1) {
      intermediate_stage_.updatePrimal(robots[omp_get_thread_num()], 
                                       primal_step_size, d[i], s[i], data_[i]);
      intermediate_stage_.updateDual(dual_step_size, data_[i]);
    }
    else {
      terminal_stage_.updatePrimal(robots[omp_get_thread_num()],  
                                   primal_step_size, d[i], s[i], data_[i]);
      terminal_stage_.updateDual(dual_step_size, data_[i]);
    }
  }
}

} // namespace robotoc