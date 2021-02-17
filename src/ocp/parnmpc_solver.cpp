#include "idocp/ocp/parnmpc_solver.hpp"

#include <stdexcept>
#include <cassert>
#include <fstream>


namespace idocp {

ParNMPCSolver::ParNMPCSolver(const Robot& robot, 
                             const std::shared_ptr<CostFunction>& cost, 
                             const std::shared_ptr<Constraints>& constraints, 
                             const double T, const int N, 
                             const int max_num_impulse, const int nthreads)
  : robots_(nthreads, robot),
    contact_sequence_(robot, N),
    parnmpc_linearizer_(N, max_num_impulse, nthreads),
    backward_correction_solver_(robot, N, max_num_impulse, nthreads),
    line_search_(robot, N, max_num_impulse, nthreads),
    parnmpc_(robot, cost, constraints, T, N, max_num_impulse),
    backward_correction_(robot, N, max_num_impulse),
    kkt_matrix_(robot, N, max_num_impulse),
    kkt_residual_(robot, N, max_num_impulse),
    s_(robot, N, max_num_impulse),
    d_(robot, N, max_num_impulse),
    N_(N),
    nthreads_(nthreads) {
  try {
    if (T <= 0) {
      throw std::out_of_range("invalid value: T must be positive!");
    }
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
    if (max_num_impulse < 0) {
      throw std::out_of_range("invalid value: max_num_impulse must be non-negative!");
    }
    if (nthreads <= 0) {
      throw std::out_of_range("invalid value: nthreads must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  for (auto& e : s_.data)    { robot.normalizeConfiguration(e.q); }
  for (auto& e : s_.impulse) { robot.normalizeConfiguration(e.q); }
  for (auto& e : s_.aux)     { robot.normalizeConfiguration(e.q); }
  for (auto& e : s_.lift)    { robot.normalizeConfiguration(e.q); }
  initConstraints();
}


ParNMPCSolver::ParNMPCSolver() {
}


ParNMPCSolver::~ParNMPCSolver() {
}


void ParNMPCSolver::initConstraints() {
  parnmpc_linearizer_.initConstraints(parnmpc_, robots_, contact_sequence_, s_);
}


void ParNMPCSolver::initBackwardCorrection(const double t) {
  parnmpc_.discretize(contact_sequence_, t);
  backward_correction_solver_.initAuxMat(parnmpc_, robots_, s_, kkt_matrix_);
}


void ParNMPCSolver::updateSolution(const double t, const Eigen::VectorXd& q, 
                                   const Eigen::VectorXd& v, 
                                   const bool use_line_search) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  parnmpc_.discretize(contact_sequence_, t);
  discretizeSolution();
  backward_correction_solver_.coarseUpdate(parnmpc_, backward_correction_, 
                                           robots_, contact_sequence_, q, v, s_, 
                                           kkt_matrix_, kkt_residual_);
  backward_correction_solver_.backwardCorrectionSerial(parnmpc_, 
                                                       backward_correction_, s_);
  backward_correction_solver_.backwardCorrectionParallel(parnmpc_, 
                                                         backward_correction_, 
                                                         robots_);
  backward_correction_solver_.forwardCorrectionSerial(parnmpc_, 
                                                      backward_correction_, 
                                                      robots_, s_);
  backward_correction_solver_.forwardCorrectionParallel(parnmpc_, 
                                                        backward_correction_, 
                                                        robots_, kkt_matrix_, 
                                                        kkt_residual_, s_, d_);
  double primal_step_size = backward_correction_solver_.maxPrimalStepSize();
  const double dual_step_size = backward_correction_solver_.maxDualStepSize();
  if (use_line_search) {
    const double max_primal_step_size = primal_step_size;
    primal_step_size = line_search_.computeStepSize(parnmpc_, robots_, 
                                                    contact_sequence_, q, v, s_,
                                                    d_, max_primal_step_size);
  }
  parnmpc_linearizer_.integrateSolution(parnmpc_, robots_, kkt_matrix_, 
                                        kkt_residual_, primal_step_size, 
                                        dual_step_size, d_, s_);
} 


const SplitSolution& ParNMPCSolver::getSolution(const int stage) const {
  assert(stage >= 0);
  assert(stage < N_);
  return s_[stage];
}


void ParNMPCSolver::getStateFeedbackGain(const int time_stage, 
                                         Eigen::MatrixXd& Kq, 
                                         Eigen::MatrixXd& Kv) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  assert(Kq.rows() == robots_[0].dimv());
  assert(Kq.cols() == robots_[0].dimv());
  assert(Kv.rows() == robots_[0].dimv());
  assert(Kv.cols() == robots_[0].dimv());
}


void ParNMPCSolver::setSolution(const std::string& name, 
                                const Eigen::VectorXd& value) {
  try {
    if (name == "q") {
      for (auto& e : s_.data)    { e.q = value; }
      for (auto& e : s_.impulse) { e.q = value; }
      for (auto& e : s_.aux)     { e.q = value; }
      for (auto& e : s_.lift)    { e.q = value; }
    }
    else if (name == "v") {
      for (auto& e : s_.data)    { e.v = value; }
      for (auto& e : s_.impulse) { e.v = value; }
      for (auto& e : s_.aux)     { e.v = value; }
      for (auto& e : s_.lift)    { e.v = value; }
    }
    else if (name == "a") {
      for (auto& e : s_.data)    { e.a  = value; }
      for (auto& e : s_.impulse) { e.dv = value; }
      for (auto& e : s_.aux)     { e.a  = value; }
      for (auto& e : s_.lift)    { e.a  = value; }
    }
    else if (name == "f") {
      for (auto& e : s_.data) { 
        for (auto& ef : e.f) { ef = value; } 
        e.set_f_stack(); 
      }
      for (auto& e : s_.impulse) { 
        for (auto& ef : e.f) { ef = value; } 
        e.set_f_stack(); 
      }
      for (auto& e : s_.aux) { 
        for (auto& ef : e.f) { ef = value; } 
        e.set_f_stack(); 
      }
      for (auto& e : s_.lift) { 
        for (auto& ef : e.f) { ef = value; } 
        e.set_f_stack(); 
      }
    }
    else if (name == "u") {
      for (auto& e : s_.data)    { e.u = value; }
      for (auto& e : s_.aux)     { e.u = value; }
      for (auto& e : s_.lift)    { e.u = value; }
    }
    else {
      throw std::invalid_argument("invalid arugment: name must be q, v, a, f, or u!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  initConstraints();
}


void ParNMPCSolver::setContactStatusUniformly(
    const ContactStatus& contact_status) {
  contact_sequence_.setContactStatusUniformly(contact_status);
  for (int i=0; i<N_; ++i) {
    s_[i].setContactStatus(contact_sequence_.contactStatus(0));
  }
}


void ParNMPCSolver::pushBackContactStatus(const ContactStatus& contact_status, 
                                          const double switching_time, 
                                          const double t) {
  const ContactStatus& last_contact_status 
      = contact_sequence_.contactStatus(contact_sequence_.numContactPhases()-1);
  DiscreteEvent discrete_event(last_contact_status, contact_status);
  discrete_event.eventTime = switching_time;
  contact_sequence_.pushBackDiscreteEvent(discrete_event);
  parnmpc_.discretize(contact_sequence_, t);
  if (discrete_event.existImpulse()) {
    const int last_impulse_index = contact_sequence_.numImpulseEvents() - 1;
    const int time_stage_before_last_impulse 
        = parnmpc_.discrete().timeStageBeforeImpulse(last_impulse_index);
    if (time_stage_before_last_impulse >= 0) {
      s_.impulse[last_impulse_index].copyPartial(
          s_[time_stage_before_last_impulse]);
      s_.aux[last_impulse_index].copy(
          s_[time_stage_before_last_impulse]);
    }
    s_.aux[last_impulse_index].setContactStatus(
        contact_sequence_.contactStatus(
            parnmpc_.discrete().contactPhaseBeforeImpulse(last_impulse_index)));
    s_.impulse[last_impulse_index].setImpulseStatus(
        contact_sequence_.impulseStatus(last_impulse_index));
    const int last_contact_phase = contact_sequence_.numContactPhases() - 1;
    for (int i=time_stage_before_last_impulse+1; i<N_; ++i) {
      s_[i].setContactStatus(contact_sequence_.contactStatus(last_contact_phase));
    }
  }
  else {
    const int last_lift_index = contact_sequence_.numLiftEvents() - 1;
    const int time_stage_before_last_lift 
        = parnmpc_.discrete().timeStageBeforeLift(last_lift_index);
    if (time_stage_before_last_lift >= 0) {
      s_.lift[last_lift_index].copy(s_[time_stage_before_last_lift]);
    }
    s_.lift[last_lift_index].setContactStatus(
        contact_sequence_.contactStatus(
            parnmpc_.discrete().contactPhaseBeforeLift(last_lift_index)));
    const int last_contact_phase = contact_sequence_.numContactPhases() - 1;
    for (int i=time_stage_before_last_lift+1; i<N_; ++i) {
      s_[i].setContactStatus(contact_sequence_.contactStatus(last_contact_phase));
    }
  }
}


void ParNMPCSolver::setContactPoints(
    const int contact_phase, 
    const std::vector<Eigen::Vector3d>& contact_points) {
  contact_sequence_.setContactPoints(contact_phase, contact_points);
}


void ParNMPCSolver::popBackDiscreteEvent() {
  if (contact_sequence_.numDiscreteEvents() > 0) {
    const int last_discrete_event = contact_sequence_.numDiscreteEvents() - 1;
    if (contact_sequence_.isImpulseEvent(last_discrete_event)) {
      const int last_impulse_index = contact_sequence_.numImpulseEvents() - 1;
      const int time_stage_after_second_last_impulse 
          = parnmpc_.discrete().timeStageAfterImpulse(last_impulse_index-1);
      const int second_last_contact_phase = contact_sequence_.numContactPhases() - 2;
      for (int i=time_stage_after_second_last_impulse; i<N_; ++i) {
        s_[i].setContactStatus(
            contact_sequence_.contactStatus(second_last_contact_phase));
      }
    }
    else {
      const int last_lift_index = contact_sequence_.numLiftEvents() - 1;
      const int time_stage_after_second_last_lift
          = parnmpc_.discrete().timeStageAfterLift(last_lift_index-1);
      const int second_last_contact_phase = contact_sequence_.numContactPhases() - 2;
      for (int i=time_stage_after_second_last_lift; i<N_; ++i) {
        s_[i].setContactStatus(
            contact_sequence_.contactStatus(second_last_contact_phase));
      }
    }
    contact_sequence_.popBackDiscreteEvent();
  }
}


void ParNMPCSolver::popFrontDiscreteEvent() {
  if (contact_sequence_.numDiscreteEvents() > 0) {
    if (contact_sequence_.isImpulseEvent(0)) {
      const int time_stage_before_first_impulse 
          = parnmpc_.discrete().timeStageBeforeImpulse(0);
      for (int i=0; i<time_stage_before_first_impulse; ++i) {
        s_[i].setContactStatus(contact_sequence_.contactStatus(1));
      }
      for (int i=0; i<contact_sequence_.numImpulseEvents()-2; ++i) {
        s_.impulse[i].copy(s_.impulse[i+1]);
      }
      for (int i=0; i<contact_sequence_.numImpulseEvents()-2; ++i) {
        s_.aux[i].copy(s_.aux[i+1]);
      }
    }
    else {
      const int time_stage_before_first_lift 
          = parnmpc_.discrete().timeStageBeforeLift(0);
      for (int i=0; i<time_stage_before_first_lift; ++i) {
        s_[i].setContactStatus(contact_sequence_.contactStatus(1));
      }
      for (int i=0; i<contact_sequence_.numLiftEvents()-2; ++i) {
        s_.lift[i].copy(s_.lift[i+1]);
      }
    }
    contact_sequence_.popFrontDiscreteEvent();
  }
}


void ParNMPCSolver::clearLineSearchFilter() {
  line_search_.clearFilter();
}


double ParNMPCSolver::KKTError() {
  return parnmpc_linearizer_.KKTError(parnmpc_, kkt_residual_);
}


void ParNMPCSolver::computeKKTResidual(const double t, const Eigen::VectorXd& q, 
                                       const Eigen::VectorXd& v) {
  parnmpc_.discretize(contact_sequence_, t);
  discretizeSolution();
  parnmpc_linearizer_.computeKKTResidual(parnmpc_, robots_, contact_sequence_, 
                                         q, v, s_, kkt_matrix_, kkt_residual_);
}


bool ParNMPCSolver::isCurrentSolutionFeasible() {
  for (int i=0; i<N_-1; ++i) {
    const bool feasible = parnmpc_[i].isFeasible(robots_[0], s_[i]);
    if (!feasible) {
      std::cout << "INFEASIBLE at time stage " << i << std::endl;
      return false;
    }
  }
  const bool feasible = parnmpc_.terminal.isFeasible(robots_[0], s_[N_-1]);
  if (!feasible) {
    std::cout << "INFEASIBLE at the terminal stage" << std::endl;
    return false;
  }
  const int num_impulse = contact_sequence_.numImpulseEvents();
  for (int i=0; i<num_impulse; ++i) {
    const bool feasible = parnmpc_.impulse[i].isFeasible(robots_[0], s_.impulse[i]);
    if (!feasible) {
      std::cout << "INFEASIBLE at impulse " << i << std::endl;
      return false;
    }
  }
  for (int i=0; i<num_impulse; ++i) {
    const bool feasible = parnmpc_.aux[i].isFeasible(robots_[0], s_.aux[i]);
    if (!feasible) {
      std::cout << "INFEASIBLE at aux " << i << std::endl;
      return false;
    }
  }
  const int num_lift = contact_sequence_.numLiftEvents();
  for (int i=0; i<num_lift; ++i) {
    const bool feasible = parnmpc_.lift[i].isFeasible(robots_[0], s_.lift[i]);
    if (!feasible) {
      std::cout << "INFEASIBLE at lift " << i << std::endl;
      return false;
    }
  }
  return true;
}


std::vector<Eigen::VectorXd> ParNMPCSolver::getSolution(
    const std::string& name) const {
  std::vector<Eigen::VectorXd> sol;
  if (name == "q") {
    for (int i=0; i<N_; ++i) {
      sol.push_back(s_[i].q);
    }
  }
  if (name == "v") {
    for (int i=0; i<N_; ++i) {
      sol.push_back(s_[i].v);
    }
  }
  if (name == "a") {
    for (int i=0; i<N_; ++i) {
      sol.push_back(s_[i].a);
    }
  }
  if (name == "f") {
    for (int i=0; i<N_; ++i) {
      sol.push_back(s_[i].f_stack());
    }
  }
  if (name == "u") {
    for (int i=0; i<N_; ++i) {
      sol.push_back(s_[i].u);
    }
  }
  return sol;
}


void ParNMPCSolver::discretizeSolution() {
  for (int i=0; i<parnmpc_.discrete().N(); ++i) {
    s_[i].setContactStatus(
        contact_sequence_.contactStatus(parnmpc_.discrete().contactPhase(i)));
  }
  for (int i=0; i<parnmpc_.discrete().numImpulseStages(); ++i) {
    s_.impulse[i].setImpulseStatus(contact_sequence_.impulseStatus(i));
    s_.aux[i].setContactStatus(
        contact_sequence_.contactStatus(
            parnmpc_.discrete().contactPhaseBeforeImpulse(i)));
  }
  for (int i=0; i<parnmpc_.discrete().numLiftStages(); ++i) {
    s_.lift[i].setContactStatus(
        contact_sequence_.contactStatus(
            parnmpc_.discrete().contactPhaseBeforeLift(i)));
  }
}

} // namespace idocp