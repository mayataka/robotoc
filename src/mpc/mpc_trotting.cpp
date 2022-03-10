#include "robotoc/mpc/mpc_trotting.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>
#include <cmath>
#include <limits>
#include <algorithm>


namespace robotoc {

MPCTrotting::MPCTrotting(const OCP& ocp, const int nthreads)
  : foot_step_planner_(std::make_shared<TrottingFootStepPlanner>(ocp.robot())),
    contact_sequence_(std::make_shared<robotoc::ContactSequence>(
        ocp.robot(), ocp.maxNumEachDiscreteEvents())),
    ocp_solver_(ocp, contact_sequence_, SolverOptions::defaultOptions(), nthreads), 
    solver_options_(SolverOptions::defaultOptions()),
    cs_standing_(ocp.robot().createContactStatus()),
    cs_lfrh_(ocp.robot().createContactStatus()),
    cs_rflh_(ocp.robot().createContactStatus()),
    vcom_(Eigen::Vector3d::Zero()),
    step_length_(Eigen::Vector3d::Zero()),
    step_height_(0),
    swing_time_(0),
    initial_lift_time_(0),
    t_(0),
    T_(ocp.T()),
    dt_(ocp.T()/ocp.N()),
    dtm_(ocp.T()/ocp.N()),
    ts_last_(0),
    eps_(std::sqrt(std::numeric_limits<double>::epsilon())),
    N_(ocp.N()),
    current_step_(0),
    predict_step_(0) {
  try {
    if (ocp.robot().maxNumPointContacts() < 4) {
      throw std::out_of_range(
          "invalid argument: robot is not a quadrupedal robot!\n robot.maxNumPointContacts() must be larger than 4!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  cs_standing_.activateContacts({0, 1, 2, 3});
  cs_lfrh_.activateContacts({0, 3});
  cs_rflh_.activateContacts({1, 2});
}


MPCTrotting::MPCTrotting() {
}


MPCTrotting::~MPCTrotting() {
}


void MPCTrotting::setGaitPattern(const Eigen::Vector3d& vcom, 
                                 const double yaw_rate, const double swing_time,
                                 const double initial_lift_time) {
  try {
    if (swing_time <= 0) {
      throw std::out_of_range("invalid value: swing_time must be positive!");
    }
    if (initial_lift_time <= 0) {
      throw std::out_of_range("invalid value: initial_lift_time must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  vcom_ = vcom;
  step_length_ = vcom * swing_time;
  swing_time_ = swing_time;
  initial_lift_time_ = initial_lift_time;
  foot_step_planner_->setGaitPattern(step_length_, (swing_time*yaw_rate));
}


void MPCTrotting::init(const double t, const Eigen::VectorXd& q, 
                       const Eigen::VectorXd& v, 
                       const SolverOptions& solver_options) {
  try {
    if (t >= initial_lift_time_) {
      throw std::out_of_range(
          "invalid value: t must be less than" + std::to_string(initial_lift_time_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  t_ = t;
  current_step_ = 0;
  predict_step_ = 0;
  contact_sequence_->initContactSequence(cs_standing_);
  bool add_step = addStep(t);
  while (add_step) {
    add_step = addStep(t);
  }
  foot_step_planner_->init(q);
  resetContactPlacements(q);
  ocp_solver_.setSolution("q", q);
  ocp_solver_.setSolution("v", v);
  ocp_solver_.setSolverOptions(solver_options);
  ocp_solver_.solve(t, q, v, true);
  ts_last_ = initial_lift_time_;
}


void MPCTrotting::setSolverOptions(const SolverOptions& solver_options) {
  ocp_solver_.setSolverOptions(solver_options);
}


void MPCTrotting::updateSolution(const double t, const double dt,
                                 const Eigen::VectorXd& q, 
                                 const Eigen::VectorXd& v) {
  assert(dt > 0);
  const bool add_step = addStep(t);
  const auto ts = contact_sequence_->eventTimes();
  bool remove_step = false;
  if (!ts.empty()) {
    if (ts.front()+eps_ < t+dt) {
      ts_last_ = ts.front();
      ocp_solver_.extrapolateSolutionInitialPhase(t);
      contact_sequence_->pop_front();
      remove_step = true;
      ++current_step_;
    }
  }
  resetContactPlacements(q);
  ocp_solver_.solve(t, q, v, true);
}


const Eigen::VectorXd& MPCTrotting::getInitialControlInput() const {
  return ocp_solver_.getSolution(0).u;
}


double MPCTrotting::KKTError(const double t, const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v) {
  return ocp_solver_.KKTError(t, q, v);
}


double MPCTrotting::KKTError() const {
  return ocp_solver_.KKTError();
}


std::shared_ptr<TrottingFootStepPlanner> MPCTrotting::getPlanner() {
  return foot_step_planner_;
}


bool MPCTrotting::addStep(const double t) {
  if (predict_step_ == 0) {
    if (initial_lift_time_ < t+T_-dtm_) {
      contact_sequence_->push_back(cs_lfrh_, initial_lift_time_);
      ++predict_step_;
      return true;
    }
  }
  else {
    double tt = ts_last_ + swing_time_;
    const auto ts = contact_sequence_->eventTimes();
    if (!ts.empty()) {
      tt = ts.back() + swing_time_;
    }
    if (tt < t+T_-dtm_) {
      if (predict_step_%2 != 0) {
        contact_sequence_->push_back(cs_rflh_, tt);
      }
      else {
        contact_sequence_->push_back(cs_lfrh_, tt);
      }
      ++predict_step_;
      return true;
    }
  }
  return false;
}


void MPCTrotting::resetContactPlacements(const Eigen::VectorXd& q) {
  const bool success = foot_step_planner_->plan(q, contact_sequence_->contactStatus(0),
                                                contact_sequence_->numContactPhases()+1);
  for (int phase=0; phase<contact_sequence_->numContactPhases(); ++phase) {
    contact_sequence_->setContactPlacements(phase, 
                                            foot_step_planner_->contactPosition(phase+1));
  }
}

} // namespace robotoc 