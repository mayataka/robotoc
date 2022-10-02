#include "robotoc/unconstr/unconstr_direct_multiple_shooting.hpp"

#include <omp.h>
#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

UnconstrDirectMultipleShooting::UnconstrDirectMultipleShooting(const OCP& ocp, 
                                                               const int nthreads)
  : nthreads_(nthreads),
    ocp_(ocp.N, SplitUnconstrOCP(ocp.robot, ocp.cost, ocp.constraints)),
    terminal_ocp_(ocp.robot, ocp.cost, ocp.constraints),
    performance_index_(),
    max_primal_step_sizes_(Eigen::VectorXd::Zero(ocp.N+1)), 
    max_dual_step_sizes_(Eigen::VectorXd::Zero(ocp.N+1)) {
  if (nthreads <= 0) {
    throw std::out_of_range("[UnconstrDirectMultipleShooting] invalid argument: nthreads must be positive!");
  }
}


UnconstrDirectMultipleShooting::UnconstrDirectMultipleShooting() {
}


void UnconstrDirectMultipleShooting::initConstraints(
    aligned_vector<Robot>& robots, 
    const std::vector<GridInfo>& time_discretization, const Solution& s) {
  const int N = time_discretization.size() - 1;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N; ++i) {
    if (i < N) {
      ocp_[i].initConstraints(robots[omp_get_thread_num()], 
                              time_discretization[i], s[i]);
    }
    else {
      terminal_ocp_.initConstraints(robots[omp_get_thread_num()], 
                                    time_discretization[i], s[i]);
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
      ocp_[i].evalOCP(robots[omp_get_thread_num()], time_discretization[i],  
                      s[i], s[i+1], kkt_residual[i]);
    }
    else {
      terminal_ocp_.evalOCP(robots[omp_get_thread_num()], time_discretization[i],  
                            s[i], kkt_residual[i]);
    }
  }
  performance_index_.setZero();
  for (int i=0; i<N; ++i) {
    performance_index_ += ocp_[i].getEval();
  }
  performance_index_ += terminal_ocp_.getEval();
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
      ocp_[i].evalKKT(robots[omp_get_thread_num()], time_discretization[i],  
                      s[i], s[i+1], kkt_matrix[i], kkt_residual[i]);
    }
    else {
      terminal_ocp_.evalKKT(robots[omp_get_thread_num()], time_discretization[i], 
                            s[i], kkt_matrix[i], kkt_residual[i]);
    }
  }
  performance_index_.setZero();
  for (int i=0; i<N; ++i) {
    performance_index_ += ocp_[i].getEval();
  }
  performance_index_ += terminal_ocp_.getEval();
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
      ocp_[i].expandPrimalAndDual(time_discretization[i].dt, kkt_matrix[i], 
                                  kkt_residual[i], d[i]);
      max_primal_step_sizes_.coeffRef(i) = ocp_[i].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(i)   = ocp_[i].maxDualStepSize();
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
      ocp_[i].updatePrimal(robots[omp_get_thread_num()], primal_step_size, 
                           d[i], s[i]);
      ocp_[i].updateDual(dual_step_size);
    }
    else {
      terminal_ocp_.updatePrimal(robots[omp_get_thread_num()],  
                                 primal_step_size, d[i], s[i]);
      terminal_ocp_.updateDual(dual_step_size);
    }
  }
}

} // namespace robotoc