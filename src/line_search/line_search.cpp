#include "robotoc/line_search/line_search.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

LineSearch::LineSearch(const OCP& ocp, 
                       const LineSearchSettings& line_search_settings) 
  : filter_(),
    settings_(line_search_settings),
    s_trial_(ocp.N+1+ocp.reserved_num_discrete_events, SplitSolution(ocp.robot)), 
    kkt_residual_(ocp.N+1+ocp.reserved_num_discrete_events, SplitKKTResidual(ocp.robot)) {
}


LineSearch::LineSearch() 
  : filter_(),
    settings_(),
    s_trial_(), 
    kkt_residual_() {
}


double LineSearch::computeStepSize(
    DirectMultipleShooting& dms, aligned_vector<Robot>& robots,
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


void LineSearch::clearFilter() {
  filter_.clear();
}


bool LineSearch::isFilterEmpty() const {
  return filter_.isEmpty();
}


void LineSearch::computeSolutionTrial(const Robot& robot, const SplitSolution& s, 
                                      const SplitDirection& d, 
                                      const double step_size, 
                                      SplitSolution& s_trial,
                                      const bool impulse) {
  s_trial.setContactStatus(s);
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_trial.q);
  s_trial.v = s.v + step_size * d.dv();
  if (!impulse) {
    s_trial.a = s.a + step_size * d.da();
    s_trial.u = s.u + step_size * d.du;
  }
  else {
    s_trial.dv = s.dv + step_size * d.ddv();
  }
  if (s.dimf() > 0) {
    s_trial.f_stack() = s.f_stack() + step_size * d.df();
    s_trial.set_f_vector();
  }
}


void LineSearch::computeSolutionTrial(const aligned_vector<Robot>& robots, 
                                      const TimeDiscretization& time_discretization, 
                                      const Solution& s, const Direction& d, 
                                      const double step_size) {
  assert(step_size > 0);
  assert(step_size <= 1);
  const int N = time_discretization.size() - 1;
  for (int i=0; i<=N; ++i) {
    computeSolutionTrial(robots[0], s[i], d[i], step_size, s_trial_[i],
                         (time_discretization[i].type == GridType::Impulse));
  }
}


double LineSearch::lineSearchFilterMethod(
    DirectMultipleShooting& dms, aligned_vector<Robot>& robots,
    const TimeDiscretization& time_discretization,
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
    const Direction& d, const double max_primal_step_size) {
  if (filter_.isEmpty()) {
    const double cost = dms.getEval().cost;
    const double violation = dms.getEval().primal_feasibility;
    filter_.augment(cost, violation);
  }
  double primal_step_size = max_primal_step_size;
  while (primal_step_size > settings_.min_step_size) {
    computeSolutionTrial(robots, time_discretization, s, d, primal_step_size);
    dms.evalOCP(robots, time_discretization, q, v, s_trial_, kkt_residual_);
    const double cost = dms.getEval().cost;
    const double violation = dms.getEval().primal_feasibility;
    if (filter_.isAccepted(cost, violation)) {
      filter_.augment(cost, violation);
      break;
    }
    primal_step_size *= settings_.step_size_reduction_rate;
  }
  return primal_step_size;
}


double LineSearch::meritBacktrackingLineSearch(
    DirectMultipleShooting& dms, aligned_vector<Robot>& robots,
    const TimeDiscretization& time_discretization,
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
    const Direction& d, const double max_primal_step_size) {
  const double penalty_param = penaltyParam(time_discretization, s);
  const double merit_now = penalty_param * dms.getEval().cost + dms.getEval().primal_feasibility;
  computeSolutionTrial(robots, time_discretization, s, d, settings_.eps);
  dms.evalOCP(robots, time_discretization, q, v, s_trial_, kkt_residual_);
  const double merit_eps = penalty_param * dms.getEval().cost + dms.getEval().primal_feasibility;
  const double directional_derivative =  (1.0 / settings_.eps) * (merit_eps - merit_now);
  double primal_step_size = max_primal_step_size;
  while (primal_step_size > settings_.min_step_size) {
    computeSolutionTrial(robots, time_discretization, s, d, primal_step_size);
    dms.evalOCP(robots, time_discretization, q, v, s_trial_, kkt_residual_);
    const double merit_next = penalty_param * dms.getEval().cost + dms.getEval().primal_feasibility;
    const bool armijoHolds = armijoCond(merit_now, merit_next, directional_derivative, 
                                        primal_step_size, settings_.armijo_control_rate);
    if (armijoHolds) {
      break;
    }
    primal_step_size *= settings_.step_size_reduction_rate;
  }
  return primal_step_size;
}


bool LineSearch::armijoCond(const double merit_now, const double merit_next, 
                            const double dd, const double step_size, 
                            const double armijo_control_rate) const {
  const double diff = armijo_control_rate * step_size * dd + merit_now - merit_next;
  return ((diff <= 0) ? true : false);
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