#include "idocp/line_search/line_search.hpp"

#include <stdexcept>
#include <cassert>


namespace idocp {

LineSearch::LineSearch(const Robot& robot, const int N, 
                       const int max_num_impulse, const int nthreads, 
                       const double step_size_reduction_rate,
                       const double min_step_size) 
  : filter_(),
    max_num_impulse_(max_num_impulse), 
    nthreads_(nthreads),
    step_size_reduction_rate_(step_size_reduction_rate), 
    min_step_size_(min_step_size),
    costs_(Eigen::VectorXd::Zero(N+1)), 
    costs_impulse_(Eigen::VectorXd::Zero(max_num_impulse)), 
    costs_aux_(Eigen::VectorXd::Zero(max_num_impulse)), 
    costs_lift_(Eigen::VectorXd::Zero(max_num_impulse)), 
    violations_(Eigen::VectorXd::Zero(N)), 
    violations_impulse_(Eigen::VectorXd::Zero(max_num_impulse)), 
    violations_aux_(Eigen::VectorXd::Zero(max_num_impulse)), 
    violations_lift_(Eigen::VectorXd::Zero(max_num_impulse)),
    s_trial_(robot, N, max_num_impulse), 
    kkt_residual_(robot, N, max_num_impulse) {
}


LineSearch::LineSearch() 
  : filter_(),
    max_num_impulse_(0), 
    nthreads_(0),
    step_size_reduction_rate_(0), 
    min_step_size_(0),
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


double LineSearch::computeStepSize(OCP& ocp, aligned_vector<Robot>& robots,
                                   const ContactSequence& contact_sequence, 
                                   const Eigen::VectorXd& q, 
                                   const Eigen::VectorXd& v, 
                                   const Solution& s, const Direction& d, 
                                   const double max_primal_step_size) {
  assert(max_primal_step_size > 0);
  assert(max_primal_step_size <= 1);
  if (filter_.isEmpty()) {
    computeCostAndViolation(ocp, robots, contact_sequence, q, v, s);
    filter_.augment(totalCosts(), totalViolations());
  }
  double primal_step_size = max_primal_step_size;
  while (primal_step_size > min_step_size_) {
    computeSolutionTrial(ocp, robots, s, d, primal_step_size);
    computeCostAndViolation(ocp, robots, contact_sequence, q, v, s_trial_,
                            primal_step_size);
    const double total_costs = totalCosts();
    const double total_violations = totalViolations();
    if (filter_.isAccepted(total_costs, total_violations)) {
      filter_.augment(total_costs, total_violations);
      break;
    }
    primal_step_size *= step_size_reduction_rate_;
  }
  if (primal_step_size > min_step_size_) {
    return primal_step_size;
  }
  else {
    return min_step_size_;
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
    const ContactSequence& contact_sequence, const Eigen::VectorXd& q, 
    const Eigen::VectorXd& v, const Solution& s, 
    const double primal_step_size) {
  assert(robots.size() == nthreads_);
  assert(q.size() == robots[0].dimq());
  assert(v.size() == robots[0].dimv());
  const int N = ocp.discrete().N();
  const int N_impulse = ocp.discrete().N_impulse();
  const int N_lift = ocp.discrete().N_lift();
  const int N_all = N + 1 + 2*N_impulse + N_lift;
  clearCosts();
  clearViolations();
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N) {
      costs_.coeffRef(i) = ocp[i].stageCost(robots[omp_get_thread_num()], 
                                            ocp.discrete().t(i), 
                                            ocp.discrete().dt(i), s[i], 
                                            primal_step_size);
      if (ocp.discrete().isTimeStageBeforeImpulse(i)) {
        violations_.coeffRef(i) = ocp[i].constraintViolation(
            robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(ocp.discrete().contactPhase(i)),
            ocp.discrete().t(i), ocp.discrete().dt(i), s[i], 
            s.impulse[ocp.discrete().impulseIndexAfterTimeStage(i)].q,
            s.impulse[ocp.discrete().impulseIndexAfterTimeStage(i)].v,
            kkt_residual_[i]);
      }
      else if (ocp.discrete().isTimeStageBeforeLift(i)) {
        violations_.coeffRef(i) = ocp[i].constraintViolation(
            robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(ocp.discrete().contactPhase(i)),
            ocp.discrete().t(i), ocp.discrete().dt(i), s[i], 
            s.lift[ocp.discrete().liftIndexAfterTimeStage(i)].q,
            s.lift[ocp.discrete().liftIndexAfterTimeStage(i)].v,
            kkt_residual_[i]);
      }
      if (ocp.discrete().isTimeStageBeforeImpulse(i+1)) {
        const int impulse_index  
            = ocp.discrete().impulseIndexAfterTimeStage(i+1);
        violations_.coeffRef(i) = ocp[i].constraintViolation(
            robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(ocp.discrete().contactPhase(i)),
            ocp.discrete().t(i), ocp.discrete().dt(i), s[i], s[i+1].q, s[i+1].v, 
            kkt_residual_[i], contact_sequence.impulseStatus(impulse_index), 
            ocp.discrete().dt(i+1), kkt_residual_.switching[impulse_index]);
      }
      else {
        violations_.coeffRef(i) = ocp[i].constraintViolation(
            robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(ocp.discrete().contactPhase(i)),
            ocp.discrete().t(i), ocp.discrete().dt(i), s[i], s[i+1].q, s[i+1].v, 
            kkt_residual_[i]);
      }
    }
    else if (i == N) {
      costs_.coeffRef(i) = ocp.terminal.terminalCost(robots[omp_get_thread_num()], 
                                                     ocp.discrete().t(i), s[i]);
    }
    else if (i < N+1+N_impulse) {
      const int impulse_index = i - (N+1);
      const int time_stage_before_impulse 
          = ocp.discrete().timeStageBeforeImpulse(impulse_index);
      costs_impulse_.coeffRef(impulse_index) 
          = ocp.impulse[impulse_index].stageCost(
              robots[omp_get_thread_num()], 
              ocp.discrete().t_impulse(impulse_index), s.impulse[impulse_index],
              primal_step_size);
      violations_impulse_.coeffRef(impulse_index) 
          = ocp.impulse[impulse_index].constraintViolation(
              robots[omp_get_thread_num()], 
              contact_sequence.impulseStatus(impulse_index), 
              ocp.discrete().t_impulse(impulse_index), s.impulse[impulse_index],
              s.aux[impulse_index].q, s.aux[impulse_index].v,
              kkt_residual_.impulse[impulse_index]);
    }
    else if (i < N+1+2*N_impulse) {
      const int impulse_index  = i - (N+1+N_impulse);
      const int time_stage_after_impulse 
          = ocp.discrete().timeStageAfterImpulse(impulse_index);
      costs_aux_.coeffRef(impulse_index) = ocp.aux[impulse_index].stageCost(
              robots[omp_get_thread_num()], 
              ocp.discrete().t_impulse(impulse_index),
              ocp.discrete().dt_aux(impulse_index), s.aux[impulse_index],
              primal_step_size);
      violations_aux_.coeffRef(impulse_index) 
          = ocp.aux[impulse_index].constraintViolation(
              robots[omp_get_thread_num()], 
              contact_sequence.contactStatus(
                  ocp.discrete().contactPhaseAfterImpulse(impulse_index)), 
              ocp.discrete().t_impulse(impulse_index), 
              ocp.discrete().dt_aux(impulse_index), s.aux[impulse_index],
              s[time_stage_after_impulse].q, s[time_stage_after_impulse].v,
              kkt_residual_.aux[impulse_index]);
    }
    else {
      const int lift_index = i - (N+1+2*N_impulse);
      const int time_stage_after_lift
          = ocp.discrete().timeStageAfterLift(lift_index);
      costs_lift_.coeffRef(lift_index) = ocp.lift[lift_index].stageCost(
              robots[omp_get_thread_num()], ocp.discrete().t_lift(lift_index), 
              ocp.discrete().dt_lift(lift_index), s.lift[lift_index], 
              primal_step_size);
      if (ocp.discrete().isTimeStageBeforeImpulse(time_stage_after_lift)) {
        const int impulse_index
            = ocp.discrete().impulseIndexAfterTimeStage(time_stage_after_lift);
        violations_lift_.coeffRef(lift_index) 
            = ocp.lift[lift_index].constraintViolation(
                robots[omp_get_thread_num()], 
                contact_sequence.contactStatus(
                    ocp.discrete().contactPhaseAfterLift(lift_index)), 
                ocp.discrete().t_lift(lift_index), 
                ocp.discrete().dt_lift(lift_index), s.lift[lift_index],
                s[time_stage_after_lift].q, s[time_stage_after_lift].v,
                kkt_residual_.lift[lift_index], 
                contact_sequence.impulseStatus(impulse_index), 
                ocp.discrete().dt(time_stage_after_lift),
                kkt_residual_.switching[impulse_index]);
      }
      else {
        violations_lift_.coeffRef(lift_index) 
            = ocp.lift[lift_index].constraintViolation(
                robots[omp_get_thread_num()], 
                contact_sequence.contactStatus(
                    ocp.discrete().contactPhaseAfterLift(lift_index)), 
                ocp.discrete().t_lift(lift_index), 
                ocp.discrete().dt_lift(lift_index), s.lift[lift_index],
                s[time_stage_after_lift].q, s[time_stage_after_lift].v,
                kkt_residual_.lift[lift_index]);
      }
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
  const int N = ocp.discrete().N();
  const int N_impulse = ocp.discrete().N_impulse();
  const int N_lift = ocp.discrete().N_lift();
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
                           s_trial_.aux[impulse_index]);
    }
    else {
      const int lift_index = i - (N+1+2*N_impulse);
      computeSolutionTrial(robots[omp_get_thread_num()], s.lift[lift_index], 
                           d.lift[lift_index], step_size, 
                           s_trial_.lift[lift_index]);
    }
  }
}

} // namespace idocp