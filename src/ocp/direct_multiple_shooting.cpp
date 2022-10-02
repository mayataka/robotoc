#include "robotoc/ocp/direct_multiple_shooting.hpp"

#include <omp.h>
#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc{

DirectMultipleShooting::DirectMultipleShooting(const OCPDef& ocp, const int nthreads)
  : ocp_data_(),
    intermediate_stage_(ocp.cost, ocp.constraints, ocp.contact_sequence),
    impact_stage_(ocp.cost, ocp.constraints, ocp.contact_sequence),
    terminal_stage_(ocp.cost, ocp.constraints, ocp.contact_sequence),
    performance_index_(),
    max_primal_step_sizes_(Eigen::VectorXd::Ones(ocp.N+1+3*ocp.num_reserved_discrete_events)), 
    max_dual_step_sizes_(Eigen::VectorXd::Ones(ocp.N+1+3*ocp.num_reserved_discrete_events)),
    nthreads_(nthreads) {
  ocp_data_.resize(ocp.N+1+3*ocp.num_reserved_discrete_events);
  for (int i=0; i<=ocp.N+3*ocp.num_reserved_discrete_events; ++i) {
    ocp_data_[i] = intermediate_stage_.createData(ocp.robot);
  }
}


DirectMultipleShooting::DirectMultipleShooting()
  : ocp_data_(),
    intermediate_stage_(),
    impact_stage_(),
    terminal_stage_(),
    performance_index_(),
    max_primal_step_sizes_(), 
    max_dual_step_sizes_(),
    nthreads_(0) {
}


void DirectMultipleShooting::initConstraints(
    aligned_vector<Robot>& robots, const TimeDiscretization& time_discretization, 
    const Solution& s) {
  const int N = time_discretization.size() - 1;
  while (ocp_data_.size() < N+1) {
    ocp_data_.push_back(ocp_data_.back());
  }
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N; ++i) {
    const auto& grid = time_discretization[i];
    if (grid.type == GridType::Terminal) {
      terminal_stage_.initConstraints(robots[omp_get_thread_num()], 
                                      grid, s[i], ocp_data_[i]);
    }
    else if (grid.type == GridType::Impulse) {
      impact_stage_.initConstraints(robots[omp_get_thread_num()], 
                                    grid, s[i], ocp_data_[i]);
    }
    else {
      intermediate_stage_.initConstraints(robots[omp_get_thread_num()], 
                                          grid, s[i], ocp_data_[i]);
    }
  }
}


bool DirectMultipleShooting::isFeasible(
    aligned_vector<Robot>& robots, const TimeDiscretization& time_discretization, 
    const Solution& s) {
  const int N = time_discretization.size() - 1;
  assert(ocp_data_.size() >= N+1);
  std::vector<bool> is_feasible(N+1, true);
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N; ++i) {
    const auto& grid = time_discretization[i];
    if (grid.type == GridType::Terminal) {
      is_feasible[i] = terminal_stage_.isFeasible(robots[omp_get_thread_num()], 
                                                  grid, s[i], ocp_data_[i]);
    }
    else if (grid.type == GridType::Impulse) {
      is_feasible[i] = impact_stage_.isFeasible(robots[omp_get_thread_num()], 
                                                grid, s[i], ocp_data_[i]);
    }
    else {
      is_feasible[i] = intermediate_stage_.isFeasible(robots[omp_get_thread_num()], 
                                                      grid, s[i], ocp_data_[i]);
    }
  }
  for (const auto e : is_feasible) {
    if (!e) return false;
  }
  return true;
}


void DirectMultipleShooting::evalOCP(
    aligned_vector<Robot>& robots, const TimeDiscretization& time_discretization, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
    KKTResidual& kkt_residual) {
  const int N = time_discretization.size() - 1;
  assert(ocp_data_.size() >= N+1);
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N; ++i) {
    const auto& grid = time_discretization[i];
    if (grid.type == GridType::Terminal) {
      terminal_stage_.evalOCP(robots[omp_get_thread_num()], grid, s[i],  
                              ocp_data_[i], kkt_residual[i]);
    }
    else if (grid.type == GridType::Impulse) {
      impact_stage_.evalOCP(robots[omp_get_thread_num()], grid, s[i], s[i+1], 
                            ocp_data_[i], kkt_residual[i]);
    }
    else {
      intermediate_stage_.evalOCP(robots[omp_get_thread_num()], grid, s[i], s[i+1], 
                                  ocp_data_[i], kkt_residual[i]);
    }
  }
  performance_index_.setZero();
  for (int i=0; i<=N; ++i) {
    performance_index_ += ocp_data_[i].performance_index;
  }
}


void DirectMultipleShooting::evalKKT(
    aligned_vector<Robot>& robots, const TimeDiscretization& time_discretization, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  const int N = time_discretization.size() - 1;
  assert(ocp_data_.size() >= N+1);
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N; ++i) {
    const auto& grid = time_discretization[i];
    if (grid.type == GridType::Terminal) {
      terminal_stage_.evalKKT(robots[omp_get_thread_num()], grid, s[i-1].q, s[i], 
                              ocp_data_[i], kkt_matrix[i], kkt_residual[i]);
    }
    else if (grid.type == GridType::Impulse) {
      impact_stage_.evalKKT(robots[omp_get_thread_num()], grid, s[i-1].q, s[i], s[i+1],
                            ocp_data_[i], kkt_matrix[i], kkt_residual[i]);
    }
    else if (i == 0) {
      intermediate_stage_.evalKKT(robots[omp_get_thread_num()], grid, q, s[i], s[i+1], 
                                  ocp_data_[i], kkt_matrix[i], kkt_residual[i]);
    }
    else {
      intermediate_stage_.evalKKT(robots[omp_get_thread_num()], grid, s[i-1].q, s[i], s[i+1],
                                  ocp_data_[i], kkt_matrix[i], kkt_residual[i]);
    }
  }
  performance_index_.setZero();
  for (int i=0; i<=N; ++i) {
    performance_index_ += ocp_data_[i].performance_index;
  }
}


void DirectMultipleShooting::computeInitialStateDirection(
    const Robot& robot,  const Eigen::VectorXd& q0, const Eigen::VectorXd& v0, 
    const Solution& s, Direction& d) const {
  ::robotoc::computeInitialStateDirection(robot, q0, v0, s[0], ocp_data_[0], d[0]);
}


const PerformanceIndex& DirectMultipleShooting::getEval() const {
  return performance_index_;
}


void DirectMultipleShooting::computeStepSizes(
    const TimeDiscretization& time_discretization, Direction& d) {
  const int N = time_discretization.size() - 1;
  assert(ocp_data_.size() >= N+1);
  if (max_primal_step_sizes_.size() < N+1) {
    max_primal_step_sizes_.resize(N+1);
    max_dual_step_sizes_.resize(N+1);
  }
  max_primal_step_sizes_.fill(1.0);
  max_dual_step_sizes_.fill(1.0);
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N; ++i) {
    const auto& grid = time_discretization[i];
    if (grid.type == GridType::Terminal) {
      terminal_stage_.expandPrimal(grid, ocp_data_[i], d[i]);
      max_primal_step_sizes_.coeffRef(i) = terminal_stage_.maxPrimalStepSize(ocp_data_[i]);
      max_dual_step_sizes_.coeffRef(i) = terminal_stage_.maxDualStepSize(ocp_data_[i]);
    }
    else if (grid.type == GridType::Impulse) {
      impact_stage_.expandPrimal(grid, ocp_data_[i], d[i]);
      max_primal_step_sizes_.coeffRef(i) = impact_stage_.maxPrimalStepSize(ocp_data_[i]);
      max_dual_step_sizes_.coeffRef(i) = impact_stage_.maxDualStepSize(ocp_data_[i]);
    }
    else {
      intermediate_stage_.expandPrimal(grid, ocp_data_[i], d[i]);
      max_primal_step_sizes_.coeffRef(i) = intermediate_stage_.maxPrimalStepSize(ocp_data_[i]);
      max_dual_step_sizes_.coeffRef(i) = intermediate_stage_.maxDualStepSize(ocp_data_[i]);
    }
  }
}


double DirectMultipleShooting::maxPrimalStepSize() const {
  return max_primal_step_sizes_.minCoeff();
}


double DirectMultipleShooting::maxDualStepSize() const {
  return max_dual_step_sizes_.minCoeff();
}


void DirectMultipleShooting::integrateSolution(
    const aligned_vector<Robot>& robots, 
    const TimeDiscretization& time_discretization, 
    const double primal_step_size, const double dual_step_size, 
    const KKTMatrix& kkt_matrix, Direction& d, Solution& s) {
  const int N = time_discretization.size() - 1;
  assert(ocp_data_.size() >= N+1);
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N; ++i) {
    const auto& grid = time_discretization[i];
    if (grid.type == GridType::Terminal) {
      terminal_stage_.expandDual(grid, ocp_data_[i], d[i]);
      terminal_stage_.updatePrimal(robots[omp_get_thread_num()], 
                                   primal_step_size, d[i], s[i], ocp_data_[i]);
      terminal_stage_.updateDual(dual_step_size, ocp_data_[i]);
    }
    else if (grid.type == GridType::Impulse) {
      impact_stage_.expandDual(grid, ocp_data_[i], d[i+1], d[i]);
      impact_stage_.updatePrimal(robots[omp_get_thread_num()], 
                                 primal_step_size, d[i], s[i], ocp_data_[i]);
      impact_stage_.updateDual(dual_step_size, ocp_data_[i]);
    }
    else {
      intermediate_stage_.expandDual(grid, ocp_data_[i], d[i+1], d[i]);
      intermediate_stage_.updatePrimal(robots[omp_get_thread_num()], 
                                       primal_step_size, d[i], s[i], ocp_data_[i]);
      intermediate_stage_.updateDual(dual_step_size, ocp_data_[i]);
    }
  }
}

} // namespace robotoc