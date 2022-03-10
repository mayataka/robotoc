#include "robotoc/mpc/mpc_jumping.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>
#include <cmath>
#include <limits>
#include <algorithm>


namespace robotoc {

MPCJumping::MPCJumping(const OCP& ocp, const int nthreads)
  : foot_step_planner_(),
    contact_sequence_(std::make_shared<robotoc::ContactSequence>(
        ocp.robot(), ocp.maxNumEachDiscreteEvents())),
    sto_constraints_(ocp.sto_constraints()),
    ocp_solver_(ocp, contact_sequence_, SolverOptions::defaultOptions(), nthreads), 
    solver_options_(SolverOptions::defaultOptions()),
    cs_ground_(ocp.robot().createContactStatus()),
    cs_flying_(ocp.robot().createContactStatus()),
    s_(),
    flying_time_(0),
    min_flying_time_(0),
    ground_time_(0),
    min_ground_time_(0),
    T_(ocp.T()),
    dt_(ocp.T()/ocp.N()),
    dtm_(1.5*(ocp.T()/ocp.N())),
    t_mpc_start_(0),
    eps_(std::sqrt(std::numeric_limits<double>::epsilon())),
    N_(ocp.N()),
    current_step_(0) {
  for (int i=0; i<cs_ground_.maxNumContacts(); ++i) {
    cs_ground_.activateContact(i);
  }
}


MPCJumping::MPCJumping() {
}


MPCJumping::~MPCJumping() {
}


void MPCJumping::setJumpPattern(
    const std::shared_ptr<FootStepPlannerBase>& foot_step_planner, 
    const double flying_time, const double min_flying_time, 
    const double ground_time, const double min_ground_time) {
  try {
    if (flying_time <= 0) {
      throw std::out_of_range("invalid value: flying_time must be positive!");
    }
    if (min_flying_time <= 0) {
      throw std::out_of_range("invalid value: min_flying_time must be positive!");
    }
    if (ground_time <= 0) {
      throw std::out_of_range("invalid value: ground_time must be positive!");
    }
    if (min_ground_time <= 0) {
      throw std::out_of_range("invalid value: min_ground_time must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  foot_step_planner_ = foot_step_planner;
  flying_time_ = flying_time;
  min_flying_time_ = min_flying_time;
  ground_time_ = ground_time;
  min_ground_time_ = min_ground_time;
}


void MPCJumping::init(const double t, const Eigen::VectorXd& q, 
                      const Eigen::VectorXd& v, 
                      const SolverOptions& solver_options, const bool sto) {
  current_step_ = 0;
  contact_sequence_->initContactSequence(cs_ground_);
  const double t_lift_off   = t + T_ - ground_time_ - flying_time_;
  const double t_touch_down = t + T_ - ground_time_;
  contact_sequence_->push_back(cs_flying_, t_lift_off, sto);
  contact_sequence_->push_back(cs_ground_, t_touch_down, sto);
  resetMinimumDwellTimes(t, dtm_);
  foot_step_planner_->init(q);
  resetContactPlacements(q);
  ocp_solver_.setSolution("q", q);
  ocp_solver_.setSolution("v", v);
  ocp_solver_.setSolverOptions(solver_options);
  ocp_solver_.solve(t, q, v, true);
  s_ = ocp_solver_.getSolution();
  const auto ts = contact_sequence_->eventTimes();
  ground_time_ = t + T_ - ts[1];
  flying_time_ = t + T_ - ts[0] - ground_time_;
  t_mpc_start_ = t;

}


void MPCJumping::reset(const double t, const Eigen::VectorXd& q, 
                       const Eigen::VectorXd& v, 
                       const SolverOptions& solver_options, const bool sto) {
  current_step_ = 0;
  contact_sequence_->initContactSequence(cs_ground_);
  const double t_lift_off   = t + T_ - ground_time_ - flying_time_;
  const double t_touch_down = t + T_ - ground_time_;
  contact_sequence_->push_back(cs_flying_, t_lift_off, sto);
  contact_sequence_->push_back(cs_ground_, t_touch_down, sto);
  resetMinimumDwellTimes(t, dtm_);
  foot_step_planner_->init(q);
  resetContactPlacements(q);
  ocp_solver_.setSolution(s_);
  ocp_solver_.setSolution("q", q);
  ocp_solver_.setSolution("q", v);
  ocp_solver_.setSolverOptions(solver_options);
  ocp_solver_.meshRefinement(t);
  ocp_solver_.solve(t, q, v, true);
  s_ = ocp_solver_.getSolution();
  const auto ts = contact_sequence_->eventTimes();
  ground_time_ = t + T_ - ts[1];
  flying_time_ = t + T_ - ts[0] - ground_time_;
  t_mpc_start_ = t;
}


void MPCJumping::setSolverOptions(const SolverOptions& solver_options) {
  ocp_solver_.setSolverOptions(solver_options);
}


void MPCJumping::updateSolution(const double t, const double dt,
                                const Eigen::VectorXd& q, 
                                const Eigen::VectorXd& v) {
  assert(dt > 0);
  const auto ts = contact_sequence_->eventTimes();
  bool remove_step = false;
  if (!ts.empty()) {
    if (ts.front()+eps_ < t+dt) {
      ocp_solver_.extrapolateSolutionInitialPhase(t);
      contact_sequence_->pop_front();
      remove_step = true;
      ++current_step_;
    }
  }
  resetMinimumDwellTimes(t, dt);
  resetContactPlacements(q);
  ocp_solver_.solve(t, q, v, true);
}


const Eigen::VectorXd& MPCJumping::getInitialControlInput() const {
  return ocp_solver_.getSolution(0).u;
}


double MPCJumping::KKTError(const double t, const Eigen::VectorXd& q, 
                            const Eigen::VectorXd& v) {
  return ocp_solver_.KKTError(t, q, v);
}


double MPCJumping::KKTError() const {
  return ocp_solver_.KKTError();
}


void MPCJumping::resetMinimumDwellTimes(const double t, const double min_dt) {
  const int num_switches = contact_sequence_->numDiscreteEvents();
  if (num_switches > 0) {
    std::vector<double> minimum_dwell_times;
    switch (current_step_) {
      case 0:
        minimum_dwell_times.push_back(min_dt);
        minimum_dwell_times.push_back(min_flying_time_);
        break;
      case 1:
        minimum_dwell_times.push_back(min_dt);
        break;
      default:
        // if current_step_ >= 2, num_switches == 0.
        break;
    }
    minimum_dwell_times.push_back(min_ground_time_+(t-t_mpc_start_));
    sto_constraints_->setMinimumDwellTimes(minimum_dwell_times);
  }
}


void MPCJumping::resetContactPlacements(const Eigen::VectorXd& q) {
  const bool success = foot_step_planner_->plan(q, contact_sequence_->contactStatus(0),
                                                contact_sequence_->numContactPhases());
  if (current_step_ == 0) {
    contact_sequence_->setContactPlacements(0, foot_step_planner_->contactPlacement(1));
    contact_sequence_->setContactPlacements(1, foot_step_planner_->contactPlacement(2));
    contact_sequence_->setContactPlacements(2, foot_step_planner_->contactPlacement(2));
  }
  else if (current_step_ == 1) {
    contact_sequence_->setContactPlacements(0, foot_step_planner_->contactPlacement(1));
    contact_sequence_->setContactPlacements(1, foot_step_planner_->contactPlacement(1));
  }
  else {
    contact_sequence_->setContactPlacements(0, foot_step_planner_->contactPlacement(1));
  }
}

} // namespace robotoc 