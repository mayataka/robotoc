#include "robotoc/mpc/mpc_quadrupedal_walking.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>
#include <cmath>
#include <algorithm>


namespace robotoc {

MPCQuadrupedalWalking::MPCQuadrupedalWalking(const OCP& ocp, 
                                             const int nthreads)
  : robot_(ocp.robot()),
    contact_sequence_(std::make_shared<robotoc::ContactSequence>(
        ocp.robot(), ocp.maxNumEachDiscreteEvents())),
    ocp_solver_(ocp, contact_sequence_, SolverOptions::defaultOptions(), nthreads), 
    solver_options_(SolverOptions::defaultOptions()),
    cs_standing_(ocp.robot().createContactStatus()),
    cs_lf_(ocp.robot().createContactStatus()),
    cs_lh_(ocp.robot().createContactStatus()),
    cs_rf_(ocp.robot().createContactStatus()),
    cs_rh_(ocp.robot().createContactStatus()),
    contact_positions_(),
    contact_positions_curr_(),
    contact_positions_prev_(),
    vcom_cmd_(Eigen::Vector3d::Zero()),
    step_length_(Eigen::Vector3d::Zero()),
    com_(Eigen::Vector3d::Zero()),
    com_curr_(Eigen::Vector3d::Zero()),
    com_prev_(Eigen::Vector3d::Zero()),
    R_(Eigen::Matrix3d::Identity()),
    R_yaw_cmd_(Eigen::Matrix3d::Identity()),
    step_height_(0),
    swing_time_(0),
    initial_lift_time_(0),
    t_(0),
    T_(ocp.T()),
    dt_(ocp.T()/ocp.N()),
    dtm_(1.5*(ocp.T()/ocp.N())),
    ts_last_(0),
    eps_(std::sqrt(std::numeric_limits<double>::epsilon())),
    N_(ocp.N()),
    current_step_(0),
    predict_step_(0) {
  cs_standing_.activateContacts({0, 1, 2, 3});
  cs_lf_.activateContacts({1, 2, 3});
  cs_lh_.activateContacts({0, 2, 3});
  cs_rf_.activateContacts({0, 1, 3});
  cs_rh_.activateContacts({0, 1, 2});
}


MPCQuadrupedalWalking::MPCQuadrupedalWalking() {
}


MPCQuadrupedalWalking::~MPCQuadrupedalWalking() {
}


void MPCQuadrupedalWalking::setGaitPattern(const Eigen::Vector3d& vcom_cmd, 
                                           const double yaw_rate_cmd,
                                           const double swing_time,
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
  vcom_cmd_ = vcom_cmd;
  step_length_ = vcom_cmd * swing_time;
  swing_time_ = swing_time;
  initial_lift_time_ = initial_lift_time;
  const double yaw_cmd = swing_time * yaw_rate_cmd;
  R_yaw_cmd_ << std::cos(yaw_cmd), -std::sin(yaw_cmd), 0, 
                std::sin(yaw_cmd), std::cos(yaw_cmd),  0,
                0, 0, 1;
}


void MPCQuadrupedalWalking::init(const double t, const Eigen::VectorXd& q, 
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

  current_step_ = 0;
  predict_step_ = 0;
  // Init contact status
  contact_sequence_->initContactSequence(cs_standing_);
  bool add_step = addStep(t);
  while (add_step) {
    add_step = addStep(t);
  }
  resetContactPlacements(q);
  ocp_solver_.setSolution("q", q);
  ocp_solver_.setSolution("v", v);
  Eigen::Vector3d f_init;
  f_init << 0, 0, 0.25*robot_.totalWeight();
  ocp_solver_.setSolution("f", f_init);
  ocp_solver_.setSolverOptions(solver_options);
  ocp_solver_.solve(t, q, v, true);
  ts_last_ = initial_lift_time_;
}


void MPCQuadrupedalWalking::setSolverOptions(
    const SolverOptions& solver_options) {
  ocp_solver_.setSolverOptions(solver_options);
}


void MPCQuadrupedalWalking::updateSolution(const double t, const double dt,
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
      contact_positions_prev_ = contact_positions_curr_;
      com_prev_ = com_curr_;
      contact_sequence_->pop_front();
      remove_step = true;
      ++current_step_;
    }
  }
  resetContactPlacements(q);
  ocp_solver_.solve(t, q, v, true);
}


const Eigen::VectorXd& MPCQuadrupedalWalking::getInitialControlInput() const {
  return ocp_solver_.getSolution(0).u;
}


double MPCQuadrupedalWalking::KKTError(const double t, const Eigen::VectorXd& q, 
                                       const Eigen::VectorXd& v) {
  return ocp_solver_.KKTError(t, q, v);
}


double MPCQuadrupedalWalking::KKTError() const {
  return ocp_solver_.KKTError();
}


bool MPCQuadrupedalWalking::addStep(const double t) {
  if (predict_step_ == 0) {
    if (initial_lift_time_ < t+T_-dtm_) {
      contact_sequence_->push_back(cs_rh_, initial_lift_time_);
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
      if (predict_step_%4 == 1) {
        contact_sequence_->push_back(cs_rf_, tt);
      }
      else if (predict_step_%4 == 2) {
        contact_sequence_->push_back(cs_lh_, tt);
      }
      else if (predict_step_%4 == 3) {
        contact_sequence_->push_back(cs_lf_, tt);
      }
      else { // predict_step_%4 == 0
        contact_sequence_->push_back(cs_rh_, tt);
      }
      ++predict_step_;
      return true;
    }
  }
  return false;
}


void MPCQuadrupedalWalking::resetContactPlacements(const Eigen::VectorXd& q) {
  robot_.updateFrameKinematics(q);
  R_ = Eigen::Quaterniond(q.coeff(6), q.coeff(3), q.coeff(4), q.coeff(5)).toRotationMatrix();
  contact_positions_.clear();
  for (const auto frame : robot_.pointContactFrames()) {
    contact_positions_.push_back(robot_.framePosition(frame));
  }
  com_ = robot_.CoM();
  contact_positions_curr_ = contact_positions_;
  com_curr_ = com_;
  // frames = [LF, LH, RF, RH] (0, 1, 2, 3)
  if (current_step_ == 0) {
    contact_positions_prev_ = contact_positions_curr_;
    com_prev_ = com_;
  }
  else if (current_step_ == 1) {
    // retrive the initial contact position from the first step (stance legs: LF, LH, RF)
    // RH
    contact_positions_[3] = contact_positions_prev_[3];
  }
  else if (current_step_ == 2) {
    // retrive the initial contact position from the first step (stance legs: LF, LH, RH)
    // RF
    contact_positions_[2] = contact_positions_prev_[2];
  }
  else if (current_step_%4 == 1) {
    // retrive the previous contact positions from the current step (stance legs: LF, LH, RF)
    // RH
    contact_positions_[3] = contact_positions_prev_[3];
  }
  else if (current_step_%4 == 2) {
    // retrive the previous contact positions from the current step (stance legs: LF, LH, RF)
    // RF
    contact_positions_[2] = contact_positions_prev_[2];
  }
  else if (current_step_%4 == 3) {
    // retrive the previous contact positions from the current step (stance legs: LF, RF, RH)
    // LH
    contact_positions_[1] = contact_positions_prev_[1];
  }
  else {
    // retrive the previous contact positions from the current step (stance legs: LH, RF, RH)
    // LF
    contact_positions_[0] = contact_positions_prev_[0];
  }
  for (int step=current_step_; step<=predict_step_; ++step) {
    R_ = (R_yaw_cmd_ * R_).eval();
    if (step == 0) {
      // do nothing (standing)
    }
    else if (step == 1) {
      contact_positions_[3].noalias() += 0.5 * R_ * step_length_;
    }
    else if (step == 2) {
      contact_positions_[2].noalias() += 0.5 * R_ * step_length_;
    }
    else if (step%4 == 1) {
      contact_positions_[3].noalias() += R_ * step_length_;
    }
    else if (step%4 == 2) {
      contact_positions_[2].noalias() += R_ * step_length_;
    }
    else if (step%4 == 3) {
      contact_positions_[1].noalias() += R_ * step_length_;
    }
    else {
      contact_positions_[0].noalias() += R_ * step_length_;
    }
    contact_sequence_->setContactPlacements(step-current_step_, contact_positions_);
  }
}

} // namespace robotoc 