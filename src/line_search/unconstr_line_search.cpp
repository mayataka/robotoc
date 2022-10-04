#include "robotoc/line_search/unconstr_line_search.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

UnconstrLineSearch::UnconstrLineSearch(const OCP& ocp,  
                                       const double step_size_reduction_rate, 
                                       const double min_step_size) 
  : filter_(),
    step_size_reduction_rate_(step_size_reduction_rate), 
    min_step_size_(min_step_size),
    s_trial_(ocp.N+1, SplitSolution(ocp.robot)), 
    kkt_residual_(ocp.N+1, SplitKKTResidual(ocp.robot)) {
}


UnconstrLineSearch::UnconstrLineSearch() 
  : filter_(),
    step_size_reduction_rate_(0), 
    min_step_size_(0),
    s_trial_(), 
    kkt_residual_() {
}


void UnconstrLineSearch::clearFilter() {
  filter_.clear();
}


bool UnconstrLineSearch::isFilterEmpty() const {
  return filter_.isEmpty();
}


double UnconstrLineSearch::computeStepSize(
    UnconstrDirectMultipleShooting& dms, aligned_vector<Robot>& robots, 
    const std::vector<GridInfo>& time_discretization, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
    const Solution& s, const Direction& d, const double max_primal_step_size) {
  assert(max_primal_step_size > 0);
  assert(max_primal_step_size <= 1);
  // If filter is empty, augment the current solution to the filter.
  if (filter_.isEmpty()) {
    const double cost = dms.getEval().cost;
    const double violation = dms.getEval().primal_feasibility;
    filter_.augment(cost, violation);
  }
  double primal_step_size = max_primal_step_size;
  while (primal_step_size > min_step_size_) {
    computeSolutionTrial(s, d, primal_step_size);
    dms.evalOCP(robots, time_discretization, q, v, s_trial_, kkt_residual_);
    const double cost = dms.getEval().cost;
    const double violation = dms.getEval().primal_feasibility;
    if (filter_.isAccepted(cost, violation)) {
      filter_.augment(cost, violation);
      break;
    }
    primal_step_size *= step_size_reduction_rate_;
  }
  return std::max(primal_step_size, min_step_size_);
}


double UnconstrLineSearch::computeStepSize(
    UnconstrBackwardCorrection& backward_correction, 
    aligned_vector<Robot>& robots, 
    const std::vector<GridInfo>& time_discretization, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
    const Solution& s, const Direction& d, const double max_primal_step_size) {
  assert(max_primal_step_size > 0);
  assert(max_primal_step_size <= 1);
  // If filter is empty, augment the current solution to the filter.
  if (filter_.isEmpty()) {
    const double cost = backward_correction.getEval().cost;
    const double violation = backward_correction.getEval().primal_feasibility;
    filter_.augment(cost, violation);
  }
  double primal_step_size = max_primal_step_size;
  while (primal_step_size > min_step_size_) {
    computeSolutionTrial(s, d, primal_step_size);
    backward_correction.evalOCP(robots, time_discretization, q, v, s_trial_, kkt_residual_);
    const double cost = backward_correction.getEval().cost;
    const double violation = backward_correction.getEval().primal_feasibility;
    if (filter_.isAccepted(cost, violation)) {
      filter_.augment(cost, violation);
      break;
    }
    primal_step_size *= step_size_reduction_rate_;
  }
  return std::max(primal_step_size, min_step_size_);
}


void UnconstrLineSearch::computeSolutionTrial(const Solution& s, 
                                              const Direction& d, 
                                              const double step_size) {
  assert(step_size > 0);
  assert(step_size <= 1);
  const int N = s.size() - 1;
  while (s_trial_.size() <= N) {
    s_trial_.push_back(s_trial_.back());
  }
  for (int i=0; i<=N; ++i) {
    computeSolutionTrial(s[i], d[i], step_size, s_trial_[i]);
  }
}


void UnconstrLineSearch::computeSolutionTrial(const SplitSolution& s, 
                                              const SplitDirection& d, 
                                              const double step_size, 
                                              SplitSolution& s_trial) {
  s_trial.q = s.q + step_size * d.dq();
  s_trial.v = s.v + step_size * d.dv();
  s_trial.a = s.a + step_size * d.da();
  s_trial.u = s.u + step_size * d.du;
}

} // namespace robotoc