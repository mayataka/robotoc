#include "robotoc/mpc/mpc_quadrupedal_trotting.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>
#include <cmath>
#include <algorithm>


namespace robotoc {

MPCQuadrupedalTrotting::MPCQuadrupedalTrotting(const OCP& ocp, 
                                               const int nthreads)
  : robot_(ocp.robot()),
    contact_sequence_(std::make_shared<robotoc::ContactSequence>(
        ocp.robot(), ocp.maxNumEachDiscreteEvents())),
    ocp_solver_(ocp, contact_sequence_, SolverOptions::defaultOptions(), nthreads), 
    solver_options_(SolverOptions::defaultOptions()),
    cs_standing_(ocp.robot().createContactStatus()),
    cs_lfrh_(ocp.robot().createContactStatus()),
    cs_rflh_(ocp.robot().createContactStatus()),
    contact_positions_(),
    contact_positions_local_(),
    vcom_cmd_(Eigen::Vector3d::Zero()),
    step_length_(Eigen::Vector3d::Zero()),
    R_(Eigen::Matrix3d::Identity()),
    R_yaw_rate_cmd_(Eigen::Matrix3d::Identity()),
    quat_(Eigen::Quaterniond::Identity()),
    step_height_(0),
    swing_time_(0),
    initial_lift_time_(0),
    T_(ocp.T()),
    dt_(ocp.T()/ocp.N()),
    dtm_(1.5*(ocp.T()/ocp.N())),
    ts_last_(0),
    eps_(std::sqrt(std::numeric_limits<double>::epsilon())),
    N_(ocp.N()),
    current_step_(0),
    predict_step_(0) {
  cs_standing_.activateContacts({0, 1, 2, 3});
  cs_lfrh_.activateContacts({0, 3});
  cs_rflh_.activateContacts({1, 2});
}


MPCQuadrupedalTrotting::MPCQuadrupedalTrotting() {
}


MPCQuadrupedalTrotting::~MPCQuadrupedalTrotting() {
}


void MPCQuadrupedalTrotting::setGaitPattern(const Eigen::Vector3d& vcom_cmd, 
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
  R_yaw_rate_cmd_ << std::cos(yaw_cmd), -std::sin(yaw_cmd), 0, 
                     std::sin(yaw_cmd), std::cos(yaw_cmd),  0,
                     0, 0, 1;
}


void MPCQuadrupedalTrotting::init(const double t, const Eigen::VectorXd& q, 
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


void MPCQuadrupedalTrotting::setSolverOptions(
    const SolverOptions& solver_options) {
  ocp_solver_.setSolverOptions(solver_options);
}


void MPCQuadrupedalTrotting::updateSolution(const double t, const double dt,
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
  if (add_step || remove_step) {
    ocp_solver_.initConstraints(t);
  }
  ocp_solver_.solve(t, q, v, false);
}


const Eigen::VectorXd& MPCQuadrupedalTrotting::getInitialControlInput() const {
  return ocp_solver_.getSolution(0).u;
}


double MPCQuadrupedalTrotting::KKTError(const double t, const Eigen::VectorXd& q, 
                                        const Eigen::VectorXd& v) {
  return ocp_solver_.KKTError(t, q, v);
}


double MPCQuadrupedalTrotting::KKTError() const {
  return ocp_solver_.KKTError();
}


bool MPCQuadrupedalTrotting::addStep(const double t) {
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


void MPCQuadrupedalTrotting::resetContactPlacements(const Eigen::VectorXd& q) {
  robot_.updateFrameKinematics(q);
  quat_ = Eigen::Quaterniond(q.coeff(6), q.coeff(3), q.coeff(4), q.coeff(5));
  R_ = quat_.toRotationMatrix();
  contact_positions_.clear();
  contact_positions_local_.clear();
  for (const auto frame : robot_.pointContactFrames()) {
    contact_positions_.push_back(robot_.framePosition(frame));
    contact_positions_local_.push_back(
        R_.transpose() * (robot_.framePosition(frame) - q.template head<3>()));
  }
  // frames = [LF, LH, RF, RH] (0, 1, 2, 3)
  if (current_step_ == 0) {
    // do nothing (standing)
  }
  else if (current_step_ == 1) {
    // retrive the initial contact positions from the first step (stance legs: LF, RH)
    // LH
    contact_positions_local_[1] << contact_positions_local_[3].coeff(0), // x : same as RH 
                                   contact_positions_local_[0].coeff(1), // y : same as LF
                                   contact_positions_local_[3].coeff(2); // z : same as RH
    contact_positions_[1].noalias() = R_ * contact_positions_local_[1] + q.template head<3>();
    // RF
    contact_positions_local_[2] << contact_positions_local_[0].coeff(0), // x : same as LF 
                                   contact_positions_local_[3].coeff(1), // y : same as RH
                                   contact_positions_local_[0].coeff(2); // z : same as LF
    contact_positions_[2].noalias() = R_ * contact_positions_local_[2] + q.template head<3>();
  }
  else if (current_step_%2 != 0) {
    // retrive the previous contact positions from the current step (stance legs: LF, RH)
    // LH
    contact_positions_local_[1] << contact_positions_local_[3].coeff(0), // x : same as RH 
                                   contact_positions_local_[0].coeff(1), // y : same as LF
                                   contact_positions_local_[3].coeff(2); // z : same as RH
    contact_positions_local_[1].noalias() -= 0.5*step_length_; 
    contact_positions_[1].noalias() = R_ * contact_positions_local_[1] + q.template head<3>();
    // RF
    contact_positions_local_[2] << contact_positions_local_[0].coeff(0), // x : same as LF 
                                   contact_positions_local_[3].coeff(1), // y : same as RH
                                   contact_positions_local_[0].coeff(2); // z : same as LF
    contact_positions_local_[2].noalias() -= 0.5*step_length_; 
    contact_positions_[2].noalias() = R_ * contact_positions_local_[2] + q.template head<3>();
  }
  else {
    // retrive the previous contact positions from the current step (stance legs: LH, RF)
    // LF
    contact_positions_local_[0] << contact_positions_local_[2].coeff(0), // x : same as RF
                                   contact_positions_local_[1].coeff(1), // y : same as LH
                                   contact_positions_local_[2].coeff(2); // z : same as RF
    contact_positions_local_[0].noalias() -= 0.5*step_length_; 
    contact_positions_[0].noalias() = R_ * contact_positions_local_[0] + q.template head<3>();
    // RH
    contact_positions_local_[3] << contact_positions_local_[1].coeff(0), // x : same as LH 
                                   contact_positions_local_[2].coeff(1), // y : same as RF
                                   contact_positions_local_[1].coeff(2); // z : same as LH
    contact_positions_local_[3].noalias() -= 0.5*step_length_; 
    contact_positions_[3].noalias() = R_ * contact_positions_local_[3] + q.template head<3>();
  }
  for (int step=current_step_; step<=predict_step_; ++step) {
    if (step == 0) {
      // do nothing (standing)
    }
    else if (step == 1) {
      if (current_step_ == 0) {
        contact_positions_[1].noalias() += 0.5 * R_ * step_length_;
        contact_positions_[2].noalias() += 0.5 * R_ * step_length_;
      }
    }
    else if (step%2 != 0) {
      contact_positions_[1].noalias() += R_ * step_length_;
      contact_positions_[2].noalias() += R_ * step_length_;
    }
    else {
      contact_positions_[0].noalias() += R_ * step_length_;
      contact_positions_[3].noalias() += R_ * step_length_;
    }
    contact_sequence_->setContactPlacements(step-current_step_, contact_positions_);
    R_ = (R_yaw_rate_cmd_ * R_).eval();
  }
}

} // namespace robotoc 