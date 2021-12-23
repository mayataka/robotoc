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
    cs_standing_(ocp.robot().createContactStatus()),
    cs_lf_(ocp.robot().createContactStatus()),
    cs_lh_(ocp.robot().createContactStatus()),
    cs_rf_(ocp.robot().createContactStatus()),
    cs_rh_(ocp.robot().createContactStatus()),
    contact_points_(),
    step_length_(0),
    step_height_(0),
    swing_time_(0),
    t0_(0),
    T_(ocp.T()),
    dt_(ocp.T()/ocp.N()),
    dtm_(1.5*(ocp.T()/ocp.N())),
    ts_last_(0),
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


void MPCQuadrupedalWalking::setGaitPattern(const double step_length, 
                                           const double step_height,
                                           const double swing_time,
                                           const double t0) {
  try {
    if (step_length <= 0) {
      throw std::out_of_range("invalid value: step_length must be positive!");
    }
    if (step_height <= 0) {
      throw std::out_of_range("invalid value: step_height must be positive!");
    }
    if (swing_time <= 0) {
      throw std::out_of_range("invalid value: swing_time must be positive!");
    }
    if (t0 <= 0) {
      throw std::out_of_range("invalid value: t0 must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  step_length_ = step_length;
  step_height_ = step_height;
  swing_time_ = swing_time;
  t0_ = t0;
}


void MPCQuadrupedalWalking::init(const double t, const Eigen::VectorXd& q, 
                                 const Eigen::VectorXd& v, 
                                 const SolverOptions& solver_options) {
  try {
    if (t >= t0_) {
      throw std::out_of_range(
          "invalid value: t must be less than" + std::to_string(t0_) + "!");
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
  double tt = t0_;
  while (tt < t+T_-dtm_) {
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
    tt += swing_time_;
  }
  resetContactPoints(q);

  ocp_solver_.setSolution("q", q);
  ocp_solver_.setSolution("v", v);
  Eigen::Vector3d f_init;
  f_init << 0, 0, 0.25*robot_.totalWeight();
  ocp_solver_.setSolution("f", f_init);
  ocp_solver_.setSolverOptions(solver_options);
  ocp_solver_.solve(t, q, v, true);
  ts_last_ = t0_;
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
    if (ts.front() < t+dt) {
      ts_last_ = ts.front();
      ocp_solver_.extrapolateSolutionInitialPhase(t);
      contact_sequence_->pop_front();
      remove_step = true;
      ++current_step_;
    }
  }
  resetContactPoints(q);
  if (add_step || remove_step) {
    ocp_solver_.initConstraints(t);
  }
  ocp_solver_.solve(t, q, v, false);
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
    if (t0_ < t+T_-dtm_) {
      contact_sequence_->push_back(cs_rh_, t0_);
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


void MPCQuadrupedalWalking::resetContactPoints(const Eigen::VectorXd& q) {
  robot_.updateFrameKinematics(q);
  contact_points_.clear();
  for (const auto frame : robot_.pointContactFrames()) {
    contact_points_.push_back(robot_.framePosition(frame));
  }
  // frames = [LF, LH, RF, RH] (0, 1, 2, 3)
  if (current_step_ == 0) {
    // do nothing (standing)
  }
  else if (current_step_ == 1) {
    // retrive the initial contact point from the first step (stance legs: LF, LH, RF)
    // RH
    contact_points_[3] << contact_points_[1].coeff(0), // x : same as LH
                          contact_points_[2].coeff(1), // y : same as RF
                          contact_points_[1].coeff(2); // z : same as LH
  }
  else if (current_step_ == 2) {
    // retrive the initial contact point from the first step (stance legs: LF, LH, RH)
    // RF
    contact_points_[2] << contact_points_[0].coeff(0), // x : same as LF 
                          contact_points_[3].coeff(1), // y : same as RH
                          contact_points_[0].coeff(2); // z : same as LF
  }
  else if (current_step_%4 == 1) {
    // retrive the previous contact points from the current step (stance legs: LF, LH, RF)
    // RH
    contact_points_[3] << contact_points_[1].coeff(0)-0.5*step_length_, // x : same as LH-0.5*step_length_ 
                          contact_points_[2].coeff(1), // y : same as RF
                          contact_points_[1].coeff(2); // z : same as LH
  }
  else if (current_step_%4 == 2) {
    // retrive the previous contact points from the current step (stance legs: LF, LH, RF)
    // RF
    contact_points_[2] << contact_points_[0].coeff(0)-0.5*step_length_, // x : same as LF-0.5*step_length_ 
                          contact_points_[3].coeff(1), // y : same as RH
                          contact_points_[0].coeff(2); // z : same as LF
  }
  else if (current_step_%4 == 3) {
    // retrive the previous contact points from the current step (stance legs: LF, RF, RH)
    // LH
    contact_points_[1] << contact_points_[3].coeff(0)-0.5*step_length_, // x : same as RH-0.5*step_length_
                          contact_points_[0].coeff(1), // y : same as LF
                          contact_points_[3].coeff(2); // z : same as RH
  }
  else {
    // retrive the previous contact points from the current step (stance legs: LH, RF, RH)
    // LF
    contact_points_[0] << contact_points_[2].coeff(0)-0.5*step_length_, // x : same as RF-0.5*step_length_
                          contact_points_[1].coeff(1), // y : same as LH
                          contact_points_[2].coeff(2); // z : same as RF
  }
  for (int step=current_step_; step<=predict_step_; ++step) {
    if (step == 0) {
      // do nothing (standing)
    }
    else if (step == 1) {
      contact_points_[3].coeffRef(0) += 0.5 * step_length_;
    }
    else if (step == 2) {
      contact_points_[2].coeffRef(0) += 0.5 * step_length_;
    }
    else if (step%4 == 1) {
      contact_points_[3].coeffRef(0) += step_length_;
    }
    else if (step%4 == 2) {
      contact_points_[2].coeffRef(0) += step_length_;
    }
    else if (step%4 == 3) {
      contact_points_[1].coeffRef(0) += step_length_;
    }
    else {
      contact_points_[0].coeffRef(0) += step_length_;
    }
    contact_sequence_->setContactPoints(step-current_step_, contact_points_);
  }
}

} // namespace robotoc 