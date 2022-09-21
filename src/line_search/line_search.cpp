#include "robotoc/line_search/line_search.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

LineSearch::LineSearch(const OCP& ocp, const int nthreads, 
                       const LineSearchSettings& line_search_settings) 
  : filter_(),
    settings_(line_search_settings),
    nthreads_(nthreads),
    costs_(Eigen::VectorXd::Zero(ocp.N()+1)), 
    costs_impulse_(Eigen::VectorXd::Zero(ocp.reservedNumDiscreteEvents())), 
    costs_aux_(Eigen::VectorXd::Zero(ocp.reservedNumDiscreteEvents())), 
    costs_lift_(Eigen::VectorXd::Zero(ocp.reservedNumDiscreteEvents())), 
    violations_(Eigen::VectorXd::Zero(ocp.N())), 
    violations_impulse_(Eigen::VectorXd::Zero(ocp.reservedNumDiscreteEvents())), 
    violations_aux_(Eigen::VectorXd::Zero(ocp.reservedNumDiscreteEvents())), 
    violations_lift_(Eigen::VectorXd::Zero(ocp.reservedNumDiscreteEvents())),
    s_trial_(ocp.robot(), ocp.N(), ocp.reservedNumDiscreteEvents()), 
    kkt_residual_(ocp.robot(), ocp.N(), ocp.reservedNumDiscreteEvents()) {
}


LineSearch::LineSearch() 
  : filter_(),
    settings_(),
    nthreads_(0),
    costs_(), 
    costs_impulse_(), 
    costs_aux_(), 
    costs_lift_(), 
    violations_(), 
    violations_impulse_(), 
    violations_aux_(), 
    violations_lift_(),
    s_trial_(), 
    kkt_residual_() {
}


LineSearch::~LineSearch() {
}


double LineSearch::computeStepSize(
    OCP& ocp, aligned_vector<Robot>& robots, 
    const std::shared_ptr<ContactSequence>& contact_sequence, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
    const Direction& d, const double max_primal_step_size) {
  assert(max_primal_step_size > 0);
  assert(max_primal_step_size <= 1);
  double primal_step_size = max_primal_step_size;
  reserve(ocp);
  if (settings_.line_search_method == LineSearchMethod::Filter) {
    primal_step_size = lineSearchFilterMethod(ocp, robots, contact_sequence, 
                                                q, v, s, d, primal_step_size);
  }
	else if (settings_.line_search_method == LineSearchMethod::MeritBacktracking) {
    primal_step_size = meritBacktrackingLineSearch(ocp, robots, contact_sequence, 
                                                     q, v, s, d, primal_step_size);
  }
  if (primal_step_size > settings_.min_step_size) {
    return primal_step_size;
  }
  else {
    return settings_.min_step_size;
  }
}


void LineSearch::clearFilter() {
  filter_.clear();
}


bool LineSearch::isFilterEmpty() const {
  return filter_.isEmpty();
}

void LineSearch::computeCostAndViolation(
    OCP& ocp, aligned_vector<Robot>& robots, 
    const std::shared_ptr<ContactSequence>& contact_sequence, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s) {
  assert(robots.size() == nthreads_);
  assert(q.size() == robots[0].dimq());
  assert(v.size() == robots[0].dimv());
  const int N = ocp.timeDiscretization().N();
  const int N_impulse = ocp.timeDiscretization().N_impulse();
  const int N_lift = ocp.timeDiscretization().N_lift();
  const int N_all = N + 1 + 2*N_impulse + N_lift;
  clearCosts();
  clearViolations();
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N) {
      if (ocp.timeDiscretization().isTimeStageBeforeImpulse(i)) {
        ocp[i].evalOCP(robots[omp_get_thread_num()], 
                       contact_sequence->contactStatus(ocp.timeDiscretization().contactPhase(i)),
                       ocp.timeDiscretization().gridInfo(i), s[i], 
                       s.impulse[ocp.timeDiscretization().impulseIndexAfterTimeStage(i)].q,
                       s.impulse[ocp.timeDiscretization().impulseIndexAfterTimeStage(i)].v,
                       kkt_residual_[i]);
        costs_.coeffRef(i) = ocp[i].stageCost();
        violations_.coeffRef(i) = ocp[i].constraintViolation(kkt_residual_[i]);
      }
      else if (ocp.timeDiscretization().isTimeStageBeforeLift(i)) {
        ocp[i].evalOCP(robots[omp_get_thread_num()], 
                       contact_sequence->contactStatus(ocp.timeDiscretization().contactPhase(i)),
                       ocp.timeDiscretization().gridInfo(i), s[i], 
                       s.lift[ocp.timeDiscretization().liftIndexAfterTimeStage(i)].q,
                       s.lift[ocp.timeDiscretization().liftIndexAfterTimeStage(i)].v,
                       kkt_residual_[i]);
        costs_.coeffRef(i) = ocp[i].stageCost();
        violations_.coeffRef(i) = ocp[i].constraintViolation(kkt_residual_[i]);
      }
      if (ocp.timeDiscretization().isTimeStageBeforeImpulse(i+1)) {
        const int impulse_index  
            = ocp.timeDiscretization().impulseIndexAfterTimeStage(i+1);
        ocp[i].evalOCP(robots[omp_get_thread_num()], 
                       contact_sequence->contactStatus(ocp.timeDiscretization().contactPhase(i)),
                       ocp.timeDiscretization().gridInfo(i), s[i], 
                       s[i+1].q, s[i+1].v, kkt_residual_[i], 
                       contact_sequence->impulseStatus(impulse_index), 
                       ocp.timeDiscretization().gridInfoAux(impulse_index), 
                       kkt_residual_.switching[impulse_index]);
        costs_.coeffRef(i) = ocp[i].stageCost();
        violations_.coeffRef(i) = ocp[i].constraintViolation(kkt_residual_[i], 
                                                             kkt_residual_.switching[impulse_index]);
      }
      else {
        ocp[i].evalOCP(robots[omp_get_thread_num()], 
                       contact_sequence->contactStatus(ocp.timeDiscretization().contactPhase(i)),
                       ocp.timeDiscretization().gridInfo(i), s[i], 
                       s[i+1].q, s[i+1].v, kkt_residual_[i]);
        costs_.coeffRef(i) = ocp[i].stageCost();
        violations_.coeffRef(i) = ocp[i].constraintViolation(kkt_residual_[i]);
      }
    }
    else if (i == N) {
      ocp.terminal.evalOCP(robots[omp_get_thread_num()], 
                           ocp.timeDiscretization().gridInfo(i), s[i], kkt_residual_[i]);
      costs_.coeffRef(i) = ocp.terminal.terminalCost();
    }
    else if (i < N+1+N_impulse) {
      const int impulse_index = i - (N+1);
      const int time_stage_before_impulse 
          = ocp.timeDiscretization().timeStageBeforeImpulse(impulse_index);
      ocp.impulse[impulse_index].evalOCP(robots[omp_get_thread_num()], 
                                         contact_sequence->impulseStatus(impulse_index), 
                                         ocp.timeDiscretization().gridInfoImpulse(impulse_index), 
                                         s.impulse[impulse_index], 
                                         s.aux[impulse_index].q, 
                                         s.aux[impulse_index].v,
                                         kkt_residual_.impulse[impulse_index]);
      costs_impulse_.coeffRef(impulse_index) = ocp.impulse[impulse_index].stageCost();
      violations_impulse_.coeffRef(impulse_index) 
          = ocp.impulse[impulse_index].constraintViolation(kkt_residual_.impulse[impulse_index]);
    }
    else if (i < N+1+2*N_impulse) {
      const int impulse_index  = i - (N+1+N_impulse);
      const int time_stage_after_impulse 
          = ocp.timeDiscretization().timeStageAfterImpulse(impulse_index);
      ocp.aux[impulse_index].evalOCP(robots[omp_get_thread_num()], 
                                     contact_sequence->contactStatus(
                                        ocp.timeDiscretization().contactPhaseAfterImpulse(impulse_index)), 
                                     ocp.timeDiscretization().gridInfoAux(impulse_index), 
                                     s.aux[impulse_index],
                                     s[time_stage_after_impulse].q, 
                                     s[time_stage_after_impulse].v,
                                     kkt_residual_.aux[impulse_index]);
      costs_aux_.coeffRef(impulse_index) = ocp.aux[impulse_index].stageCost();
      violations_aux_.coeffRef(impulse_index) 
          = ocp.aux[impulse_index].constraintViolation(kkt_residual_.aux[impulse_index]);
    }
    else {
      const int lift_index = i - (N+1+2*N_impulse);
      const int time_stage_after_lift
          = ocp.timeDiscretization().timeStageAfterLift(lift_index);
      ocp.lift[lift_index].evalOCP(robots[omp_get_thread_num()], 
                                   contact_sequence->contactStatus(
                                       ocp.timeDiscretization().contactPhaseAfterLift(lift_index)), 
                                   ocp.timeDiscretization().gridInfoLift(lift_index), 
                                   s.lift[lift_index],
                                   s[time_stage_after_lift].q, 
                                   s[time_stage_after_lift].v,
                                   kkt_residual_.lift[lift_index]);
      violations_lift_.coeffRef(lift_index) 
          = ocp.lift[lift_index].constraintViolation(kkt_residual_.lift[lift_index]);
      costs_lift_.coeffRef(lift_index) = ocp.lift[lift_index].stageCost();
    }
  }
}


void LineSearch::computeSolutionTrial(const OCP& ocp, 
                                      const aligned_vector<Robot>& robots, 
                                      const Solution& s, const Direction& d, 
                                      const double step_size) {
  assert(robots.size() == nthreads_);
  assert(step_size > 0);
  assert(step_size <= 1);
  const int N = ocp.timeDiscretization().N();
  const int N_impulse = ocp.timeDiscretization().N_impulse();
  const int N_lift = ocp.timeDiscretization().N_lift();
  const int N_all = N + 1 + 2*N_impulse + N_lift;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i <= N) {
      computeSolutionTrial(robots[omp_get_thread_num()], s[i], d[i], step_size, 
                           s_trial_[i]);
    }
    else if (i < N+1+N_impulse) {
      const int impulse_index = i - (N+1);
      computeSolutionTrial(robots[omp_get_thread_num()], 
                           s.impulse[impulse_index], 
                           d.impulse[impulse_index], step_size, 
                           s_trial_.impulse[impulse_index]);
    }
    else if (i < N+1+2*N_impulse) {
      const int impulse_index  = i - (N+1+N_impulse);
      computeSolutionTrial(robots[omp_get_thread_num()], s.aux[impulse_index], 
                           d.aux[impulse_index], step_size, 
                           s_trial_.aux[impulse_index], true);
    }
    else {
      const int lift_index = i - (N+1+2*N_impulse);
      computeSolutionTrial(robots[omp_get_thread_num()], s.lift[lift_index], 
                           d.lift[lift_index], step_size, 
                           s_trial_.lift[lift_index]);
    }
  }
}


double LineSearch::lineSearchFilterMethod(
    OCP& ocp, aligned_vector<Robot>& robots, 
    const std::shared_ptr<ContactSequence>& contact_sequence, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
    const Direction& d, const double initial_primal_step_size) {
  if (filter_.isEmpty()) {
    computeCostAndViolation(ocp, robots, contact_sequence, q, v, s);
    filter_.augment(totalCosts(), totalViolations());
  }
  double primal_step_size = initial_primal_step_size;
  while (primal_step_size > settings_.min_step_size) {
    computeSolutionTrial(ocp, robots, s, d, primal_step_size);
    computeCostAndViolation(ocp, robots, contact_sequence, q, v, s_trial_);
    const double total_costs = totalCosts();
    const double total_violations = totalViolations();
    if (filter_.isAccepted(total_costs, total_violations)) {
      filter_.augment(total_costs, total_violations);
      break;
    }
    primal_step_size *= settings_.step_size_reduction_rate;
  }
  return primal_step_size;
}


double LineSearch::meritBacktrackingLineSearch(
    OCP& ocp, aligned_vector<Robot>& robots, 
    const std::shared_ptr<ContactSequence>& contact_sequence, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
    const Direction& d, const double initial_primal_step_size) {
  computeCostAndViolation(ocp, robots, contact_sequence, q, v, s);
  const double penalty_param = penaltyParam(ocp, s);
  const double merit_now = merit(penalty_param);
  computeSolutionTrial(ocp, robots, s, d, settings_.eps);
  computeCostAndViolation(ocp, robots, contact_sequence, q, v, s_trial_);
  const double merit_eps = merit(penalty_param);
  const double directional_derivative =  (1.0 / settings_.eps) * (merit_eps - merit_now);
  double primal_step_size = initial_primal_step_size;
  while (primal_step_size > settings_.min_step_size) {
    computeSolutionTrial(ocp, robots, s, d, primal_step_size);
    computeCostAndViolation(ocp, robots, contact_sequence, q, v, s_trial_);
    const double merit_next = merit(penalty_param);
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


double LineSearch::penaltyParam(const OCP& ocp, const Solution& s) const {
  const int N = ocp.timeDiscretization().N();
  const int N_impulse = ocp.timeDiscretization().N_impulse();
  const int N_lift = ocp.timeDiscretization().N_lift();
  const int N_all = N + 1 + 2*N_impulse + N_lift;
  Eigen::VectorXd lagrangeMultiplierLinfNorms = Eigen::VectorXd::Zero(N_all);                                        
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i <= N) {
      lagrangeMultiplierLinfNorms.coeffRef(i) = s[i].lagrangeMultiplierLinfNorm();
    }
    else if (i < N+1+N_impulse) {
      const int impulse_index = i - (N+1);
      lagrangeMultiplierLinfNorms.coeffRef(i) = s.impulse[impulse_index].lagrangeMultiplierLinfNorm();
    }
    else if (i < N+1+2*N_impulse) {
      const int impulse_index  = i - (N+1+N_impulse);
      lagrangeMultiplierLinfNorms.coeffRef(i) = s.aux[impulse_index].lagrangeMultiplierLinfNorm();
    }
    else {
      const int lift_index = i - (N+1+2*N_impulse);
      lagrangeMultiplierLinfNorms.coeffRef(i) = s.lift[lift_index].lagrangeMultiplierLinfNorm();
    }
  }
  return lagrangeMultiplierLinfNorms.maxCoeff() * (1 + settings_.margin_rate);
}


double LineSearch::merit(const double penalty_param) const {
  const double res = (costs_.head(violations_.size()) + penalty_param * violations_).sum() + costs_[violations_.size()];
  const double res_impulse = (costs_impulse_ + penalty_param * violations_impulse_).sum();
  const double res_aux = (costs_aux_ + penalty_param * violations_aux_).sum();
  const double res_lift = (costs_lift_ + penalty_param * violations_lift_).sum();                                       
  return res + res_impulse + res_aux + res_lift;
}


void LineSearch::set(const LineSearchSettings& settings) {
  settings_ = settings;
}


void LineSearch::reserve(const OCP& ocp) {
  costs_impulse_.resize(ocp.reservedNumDiscreteEvents());
  costs_impulse_.setZero();
  costs_aux_.resize(ocp.reservedNumDiscreteEvents());
  costs_aux_.setZero();
  costs_lift_.resize(ocp.reservedNumDiscreteEvents());
  costs_lift_.setZero();
  violations_impulse_.resize(ocp.reservedNumDiscreteEvents());
  violations_impulse_.setZero();
  violations_aux_.resize(ocp.reservedNumDiscreteEvents());
  violations_aux_.setZero();
  violations_lift_.resize(ocp.reservedNumDiscreteEvents());
  violations_lift_.setZero();
  s_trial_.reserve(ocp.robot(), ocp.reservedNumDiscreteEvents());
  kkt_residual_.reserve(ocp.robot(), ocp.reservedNumDiscreteEvents());
}

} // namespace robotoc