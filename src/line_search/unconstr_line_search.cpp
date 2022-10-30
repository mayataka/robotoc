#include "robotoc/line_search/unconstr_line_search.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

UnconstrLineSearch::UnconstrLineSearch(const OCP& ocp,  
                                       const LineSearchSettings& settings) 
  : filter_(settings.filter_cost_reduction_rate, 
            settings.filter_constraint_violation_reduction_rate),
    settings_(settings),
    dms_trial_(ocp, 1),
    bc_trial_(ocp, 1),
    s_trial_(ocp.N+1, SplitSolution(ocp.robot)), 
    kkt_residual_(ocp.N+1, SplitKKTResidual(ocp.robot)) {
}


UnconstrLineSearch::UnconstrLineSearch() 
  : filter_(),
    settings_(), 
    dms_trial_(),
    bc_trial_(),
    s_trial_(), 
    kkt_residual_() {
}


void UnconstrLineSearch::clearHistory() {
  filter_.clear();
}


double UnconstrLineSearch::computeStepSize(
    const UnconstrDirectMultipleShooting& dms, aligned_vector<Robot>& robots, 
    const std::vector<GridInfo>& time_discretization, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
    const Solution& s, const Direction& d, const double max_primal_step_size) {
  assert(max_primal_step_size > 0);
  assert(max_primal_step_size <= 1);
  if (filter_.isEmpty()) {
    const double cost = dms.getEval().cost + dms.getEval().cost_barrier;
    const double violation = dms.getEval().primal_feasibility;
    filter_.augment(cost, violation);
  }
  double primal_step_size = max_primal_step_size;
  while (primal_step_size > settings_.min_step_size) {
    dms_trial_ = dms;
    s_trial_ = s;
    dms_trial_.integratePrimalSolution(robots, time_discretization, 
                                       primal_step_size, d, s_trial_);
    dms_trial_.evalOCP(robots, time_discretization, q, v, s_trial_, kkt_residual_);
    const double cost = dms_trial_.getEval().cost + dms_trial_.getEval().cost_barrier;
    const double violation = dms_trial_.getEval().primal_feasibility;
    if (filter_.isAccepted(cost, violation)) {
      filter_.augment(cost, violation);
      return primal_step_size;
    }
    primal_step_size *= settings_.step_size_reduction_rate;
  }
  return std::max(primal_step_size, settings_.min_step_size);
}


double UnconstrLineSearch::computeStepSize(
    const UnconstrBackwardCorrection& bc, aligned_vector<Robot>& robots, 
    const std::vector<GridInfo>& time_discretization, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
    const Solution& s, const Direction& d, const double max_primal_step_size) {
  assert(max_primal_step_size > 0);
  assert(max_primal_step_size <= 1);
  if (filter_.isEmpty()) {
    const double cost = bc.getEval().cost + bc.getEval().cost_barrier;
    const double violation = bc.getEval().primal_feasibility;
    filter_.augment(cost, violation);
  }
  double primal_step_size = max_primal_step_size;
  while (primal_step_size > settings_.min_step_size) {
    bc_trial_ = bc;
    s_trial_ = s;
    bc_trial_.integratePrimalSolution(robots, time_discretization, 
                                      primal_step_size, d, s_trial_);
    bc_trial_.evalOCP(robots, time_discretization, q, v, s_trial_, kkt_residual_);
    const double cost = bc_trial_.getEval().cost + bc_trial_.getEval().cost_barrier;
    const double violation = bc_trial_.getEval().primal_feasibility;
    if (filter_.isAccepted(cost, violation)) {
      filter_.augment(cost, violation);
      return primal_step_size;
    }
    primal_step_size *= settings_.step_size_reduction_rate;
  }
  return std::max(primal_step_size, settings_.min_step_size);
}

} // namespace robotoc