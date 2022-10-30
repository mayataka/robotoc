#include "robotoc/line_search/line_search.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

LineSearch::LineSearch(const OCP& ocp, 
                       const LineSearchSettings& line_search_settings) 
  : filter_(),
    settings_(line_search_settings),
    dms_trial_(ocp, 1),
    s_trial_(ocp.N+1+ocp.reserved_num_discrete_events, SplitSolution(ocp.robot)), 
    kkt_residual_(ocp.N+1+ocp.reserved_num_discrete_events, SplitKKTResidual(ocp.robot)) {
}


LineSearch::LineSearch() 
  : filter_(),
    settings_(),
    dms_trial_(),
    s_trial_(), 
    kkt_residual_() {
}


double LineSearch::computeStepSize(
    const DirectMultipleShooting& dms, aligned_vector<Robot>& robots,
    const TimeDiscretization& time_discretization,
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
    const Direction& d, const double max_primal_step_size) {
  assert(max_primal_step_size > 0);
  assert(max_primal_step_size <= 1);
  double primal_step_size = max_primal_step_size;
  resizeData(time_discretization);
  if (settings_.line_search_method == LineSearchMethod::Filter) {
    primal_step_size = lineSearchFilterMethod(dms, robots, time_discretization, 
                                              q, v, s, d, primal_step_size);
  }
	else if (settings_.line_search_method == LineSearchMethod::MeritBacktracking) {
    primal_step_size = meritBacktrackingLineSearch(dms, robots, time_discretization, 
                                                   q, v, s, d, primal_step_size);
  }
  else {
    throw std::runtime_error("[LineSearch]: Invalid LineSearchMethod");
  }
  return std::max(primal_step_size, settings_.min_step_size);
}


void LineSearch::clearHistory() {
  filter_.clear();
}


bool LineSearch::isFilterEmpty() const {
  return filter_.isEmpty();
}


double LineSearch::lineSearchFilterMethod(
    const DirectMultipleShooting& dms, aligned_vector<Robot>& robots,
    const TimeDiscretization& time_discretization,
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
    const Direction& d, const double max_primal_step_size) {
  if (filter_.isEmpty()) {
    const double cost = dms_trial_.getEval().cost;
    const double violation = dms_trial_.getEval().primal_feasibility;
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
  return primal_step_size;
}


double LineSearch::meritBacktrackingLineSearch(
    const DirectMultipleShooting& dms, aligned_vector<Robot>& robots,
    const TimeDiscretization& time_discretization,
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
    const Direction& d, const double max_primal_step_size) {
  const double penalty_param = penaltyParam(time_discretization, s);
  const double merit = penalty_param * dms.getEval().cost + dms.getEval().primal_feasibility;
  dms_trial_ = dms;
  s_trial_ = s;
  dms_trial_.integratePrimalSolution(robots, time_discretization, settings_.eps, d, s_trial_);
  dms_trial_.evalOCP(robots, time_discretization, q, v, s_trial_, kkt_residual_);
  const double merit_eps = dms_trial_.getEval().cost + dms_trial_.getEval().cost_barrier 
                            + penalty_param * dms_trial_.getEval().primal_feasibility;
  const double merit_directional_derivative = (1.0 / settings_.eps) * (merit_eps - merit);

  double primal_step_size = max_primal_step_size;
  while (primal_step_size > settings_.min_step_size) {
    dms_trial_ = dms;
    s_trial_ = s;
    dms_trial_.integratePrimalSolution(robots, time_discretization, 
                                       primal_step_size, d, s_trial_);
    dms_trial_.evalOCP(robots, time_discretization, q, v, s_trial_, kkt_residual_);
    const double merit_trial = dms_trial_.getEval().cost + dms_trial_.getEval().cost_barrier 
                                + penalty_param * dms_trial_.getEval().primal_feasibility;
    if (armijoCondition(merit, merit_trial, merit_directional_derivative, primal_step_size)) {
      return primal_step_size;
    }
    primal_step_size *= settings_.step_size_reduction_rate;
  }
  return primal_step_size;
}


bool LineSearch::armijoCondition(const double merit, const double merit_trial, 
                                 const double merit_directional_derivative, 
                                 const double step_size) const {
  const double merit_linear_model 
      = merit + settings_.armijo_control_rate * step_size * merit_directional_derivative;
  return (merit_trial < merit_linear_model);
}


double LineSearch::penaltyParam(const TimeDiscretization time_discretization, 
                                const Solution& s) const {
  const int N = time_discretization.size() - 1;
  Eigen::VectorXd lagrangeMultiplierLinfNorms = Eigen::VectorXd::Zero(N+1);                                        
  for (int i=0; i<=N; ++i) {
    lagrangeMultiplierLinfNorms.coeffRef(i) = s[i].lagrangeMultiplierLinfNorm();
  }
  return lagrangeMultiplierLinfNorms.maxCoeff() * (1 + settings_.margin_rate);
}


void LineSearch::set(const LineSearchSettings& settings) {
  settings_ = settings;
}


void LineSearch::resizeData(const TimeDiscretization& time_discretization) {
  while (s_trial_.size() < time_discretization.size()) {
    s_trial_.push_back(s_trial_.back());
  }
  while (kkt_residual_.size() < time_discretization.size()) {
    kkt_residual_.push_back(kkt_residual_.back());
  }
}

} // namespace robotoc