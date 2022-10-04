#include "robotoc/unconstr/unconstr_direct_multiple_shooting.hpp"

#include <omp.h>
#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

UnconstrDirectMultipleShooting::UnconstrDirectMultipleShooting(const OCP& ocp, 
                                                               const int nthreads)
  : nthreads_(nthreads),
    intermediate_stage_(ocp.robot, ocp.cost, ocp.constraints),
    terminal_stage_(ocp.robot, ocp.cost, ocp.constraints),
    data_(),
    performance_index_(),
    max_primal_step_sizes_(Eigen::VectorXd::Zero(ocp.N+1)), 
    max_dual_step_sizes_(Eigen::VectorXd::Zero(ocp.N+1)) {
  if (nthreads <= 0) {
    throw std::out_of_range("[UnconstrDirectMultipleShooting] invalid argument: nthreads must be positive!");
  }
  data_.resize(ocp.N+1);
  for (int i=0; i<ocp.N; ++i) {
    data_[i] = intermediate_stage_.createData(ocp.robot);
  }
  data_[ocp.N] = terminal_stage_.createData(ocp.robot);
}


UnconstrDirectMultipleShooting::UnconstrDirectMultipleShooting() {
}


void UnconstrDirectMultipleShooting::setNumThreads(const int nthreads) {
  if (nthreads <= 0) {
    throw std::out_of_range("[UnconstrDirectMultipleShooting] invalid argument: nthreads must be positive!");
  }
  nthreads_ = nthreads;
}


void UnconstrDirectMultipleShooting::initConstraints(
    aligned_vector<Robot>& robots, 
    const std::vector<GridInfo>& time_discretization, const Solution& s) {
  const int N = time_discretization.size() - 1;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N; ++i) {
    if (i < N) {
      intermediate_stage_.initConstraints(robots[omp_get_thread_num()], 
                                          time_discretization[i], s[i], data_[i]);
    }
    else {
      terminal_stage_.initConstraints(robots[omp_get_thread_num()], 
                                      time_discretization[i], s[i], data_[i]);
    }
  }
}


void UnconstrDirectMultipleShooting::evalOCP(
    aligned_vector<Robot>& robots, 
    const std::vector<GridInfo>& time_discretization, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
    KKTResidual& kkt_residual) {
  assert(robots.size() >= nthreads_);
  const int N = time_discretization.size() - 1;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N; ++i) {
    if (i < N) {
      intermediate_stage_.evalOCP(robots[omp_get_thread_num()], 
                                  time_discretization[i],  
                                  s[i], s[i+1], data_[i], kkt_residual[i]);
    }
    else {
      terminal_stage_.evalOCP(robots[omp_get_thread_num()], 
                              time_discretization[i], 
                              s[i], data_[i], kkt_residual[i]);
    }
  }
  performance_index_.setZero();
  for (int i=0; i<=N; ++i) {
    performance_index_ += data_[i].performance_index;
  }
}


void UnconstrDirectMultipleShooting::evalKKT(
    aligned_vector<Robot>& robots, 
    const std::vector<GridInfo>& time_discretization, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  assert(robots.size() >= nthreads_);
  const int N = time_discretization.size() - 1;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N; ++i) {
    if (i < N) {
      intermediate_stage_.evalKKT(robots[omp_get_thread_num()], 
                                  time_discretization[i], s[i], s[i+1], 
                                  data_[i], kkt_matrix[i], kkt_residual[i]);
    }
    else {
      terminal_stage_.evalKKT(robots[omp_get_thread_num()], 
                              time_discretization[i], s[i], 
                              data_[i], kkt_matrix[i], kkt_residual[i]);
    }
  }
  performance_index_.setZero();
  for (int i=0; i<=N; ++i) {
    performance_index_ += data_[i].performance_index;
  }
}


void UnconstrDirectMultipleShooting::computeInitialStateDirection(
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
    Direction& d) {
  d[0].dq() = q - s[0].q;
  d[0].dv() = v - s[0].v;
}


const PerformanceIndex& UnconstrDirectMultipleShooting::getEval() const {
  return performance_index_;
}


void UnconstrDirectMultipleShooting::computeStepSizes(
    const std::vector<GridInfo>& time_discretization, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual, Direction& d) {
  const int N = time_discretization.size() - 1;
  max_primal_step_sizes_.fill(1.0);
  max_dual_step_sizes_.fill(1.0);
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N; ++i) {
    if (i < N) {
      intermediate_stage_.expandPrimalAndDual(time_discretization[i].dt, 
                                              kkt_matrix[i], kkt_residual[i], 
                                              data_[i], d[i]);
      max_primal_step_sizes_.coeffRef(i) 
          = intermediate_stage_.maxPrimalStepSize(data_[i]);
      max_dual_step_sizes_.coeffRef(i)
          = intermediate_stage_.maxDualStepSize(data_[i]);
    }
  }
}


double UnconstrDirectMultipleShooting::maxPrimalStepSize() const {
  return max_primal_step_sizes_.minCoeff();
}


double UnconstrDirectMultipleShooting::maxDualStepSize() const {
  return max_dual_step_sizes_.minCoeff();
}


void UnconstrDirectMultipleShooting::integrateSolution(
    const aligned_vector<Robot>& robots,
    const std::vector<GridInfo>& time_discretization, 
    const double primal_step_size, const double dual_step_size,
    Direction& d, Solution& s) {
  assert(robots.size() >= nthreads_);
  const int N = time_discretization.size() - 1;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N; ++i) {
    if (i < N) {
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