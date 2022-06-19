#include "robotoc/mpc/trot_foot_step_planner.hpp"
#include "robotoc/utils/rotation.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

TrotFootStepPlanner::TrotFootStepPlanner(const Robot& quadruped_robot)
  : ContactPlannerBase(),
    robot_(quadruped_robot),
    raibert_heuristic_(),
    vcom_moving_window_filter_(),
    enable_raibert_heuristic_(false),
    LF_foot_id_(quadruped_robot.pointContactFrames()[0]),
    LH_foot_id_(quadruped_robot.pointContactFrames()[1]),
    RF_foot_id_(quadruped_robot.pointContactFrames()[2]),
    RH_foot_id_(quadruped_robot.pointContactFrames()[3]),
    current_step_(0),
    contact_position_ref_(),
    com_ref_(),
    R_(),
    com_to_contact_position_local_(),
    vcom_(Eigen::Vector3d::Zero()),
    vcom_cmd_(Eigen::Vector3d::Zero()),
    step_length_(Eigen::Vector3d::Zero()),
    R_yaw_(Eigen::Matrix3d::Identity()),
    R_current_(Eigen::Matrix3d::Identity()),
    yaw_rate_cmd_(0),
    enable_stance_phase_(false) {
  try {
    if (quadruped_robot.maxNumPointContacts() < 4) {
      throw std::out_of_range(
          "invalid argument: robot is not a quadrupedal robot!\n robot.maxNumPointContacts() must be larger than 4!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


TrotFootStepPlanner::TrotFootStepPlanner() {
}


TrotFootStepPlanner::~TrotFootStepPlanner() {
}


void TrotFootStepPlanner::setGaitPattern(const Eigen::Vector3d& step_length, 
                                         const double step_yaw, 
                                         const bool enable_stance_phase) {
  step_length_ = step_length;
  R_yaw_<< std::cos(step_yaw), -std::sin(step_yaw), 0, 
           std::sin(step_yaw), std::cos(step_yaw),  0,
           0, 0, 1;
  enable_stance_phase_ = enable_stance_phase;
  enable_raibert_heuristic_ = false;
}


void TrotFootStepPlanner::setRaibertGaitPattern(const Eigen::Vector3d& vcom_cmd,
                                                const double yaw_rate_cmd, 
                                                const double swing_time,
                                                const double stance_time,
                                                const double gain) {
  try {
    if (swing_time <= 0.0) {
      throw std::out_of_range("invalid argument: swing_time must be positive!");
    }
    if (stance_time < 0.0) {
      throw std::out_of_range("invalid argument: stance_time must be non-negative!");
    }
    if (gain <= 0.0) {
      throw std::out_of_range("invalid argument: gain must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  const double period = 2.0 * (swing_time + stance_time);
  raibert_heuristic_.setParameters(period, gain);
  vcom_moving_window_filter_.setParameters(period, 0.1*period);
  vcom_cmd_ = vcom_cmd;
  const double step_yaw = yaw_rate_cmd * (swing_time + stance_time);
  R_yaw_<< std::cos(step_yaw), -std::sin(step_yaw), 0, 
           std::sin(step_yaw),  std::cos(step_yaw), 0,
           0, 0, 1;
  yaw_rate_cmd_ = yaw_rate_cmd;
  enable_stance_phase_ = (stance_time > 0.0);
  enable_raibert_heuristic_ = true;
}


void TrotFootStepPlanner::init(const Eigen::VectorXd& q) {
  Eigen::Matrix3d R = rotation::toRotationMatrix(q.template segment<4>(3));
  rotation::projectRotationMatrix(R, rotation::ProjectionAxis::Z);
  robot_.updateFrameKinematics(q);
  com_to_contact_position_local_ = { R.transpose() * (robot_.framePosition(LF_foot_id_)-robot_.CoM()), 
                                     R.transpose() * (robot_.framePosition(LH_foot_id_)-robot_.CoM()),
                                     R.transpose() * (robot_.framePosition(RF_foot_id_)-robot_.CoM()),
                                     R.transpose() * (robot_.framePosition(RH_foot_id_)-robot_.CoM()) };
  contact_position_ref_.clear();
  com_ref_.clear(),
  com_ref_.push_back(robot_.CoM());
  R_.clear();
  R_.push_back(R);
  vcom_moving_window_filter_.clear();
  current_step_ = 0;
}


bool TrotFootStepPlanner::plan(const double t, const Eigen::VectorXd& q,
                               const Eigen::VectorXd& v,
                               const ContactStatus& contact_status,
                               const int planning_steps) {
  assert(planning_steps >= 0);
  if (enable_raibert_heuristic_) {
    vcom_ = v.template head<3>();
    vcom_moving_window_filter_.push_back(t, vcom_.template head<2>());
    const Eigen::Vector2d& vcom_avg = vcom_moving_window_filter_.average();
    raibert_heuristic_.planStepLength(vcom_avg, vcom_cmd_.template head<2>(), 
                                      yaw_rate_cmd_);
    step_length_ = raibert_heuristic_.stepLength();
  }
  robot_.updateFrameKinematics(q);
  std::vector<Eigen::Vector3d> contact_position;
  for (const auto frame : robot_.pointContactFrames()) {
    contact_position.push_back(robot_.framePosition(frame));
  }
  Eigen::Vector3d com = Eigen::Vector3d::Zero();
  Eigen::Matrix3d R = R_.front();
  if (contact_status.isContactActive(0) && contact_status.isContactActive(1) 
      && contact_status.isContactActive(2) && contact_status.isContactActive(3)) {
    if (enable_stance_phase_) {
      if (current_step_%2 != 0) {
        ++current_step_;
      }
    }
    else {
      current_step_ = 0;
    }
    for (int i=0; i<4; ++i) {
      com.noalias() += contact_position[i];
      com.noalias() -= R * com_to_contact_position_local_[i];
    }
    com.array() /= 4.0;
  }
  else if (contact_status.isContactActive(0) && contact_status.isContactActive(3)) {
    if (enable_stance_phase_) {
      if (current_step_%4 != 1) {
        ++current_step_;
        R = (R_yaw_ * R).eval();
      }
    }
    else {
      if (current_step_%2 != 1) {
        ++current_step_;
        R = (R_yaw_ * R).eval();
      }
    }
    com.noalias() += contact_position[0];
    com.noalias() -= R * com_to_contact_position_local_[0];
    com.noalias() += contact_position[3];
    com.noalias() -= R * com_to_contact_position_local_[3];
    com.array() /= 2.0;
    contact_position[1].noalias() = com + R * (com_to_contact_position_local_[1] - 0.5 * step_length_);
    contact_position[2].noalias() = com + R * (com_to_contact_position_local_[2] - 0.5 * step_length_);
  }
  else if (contact_status.isContactActive(1) && contact_status.isContactActive(2)) {
    if (enable_stance_phase_) {
      if (current_step_%4 != 3) {
        ++current_step_;
        R = (R_yaw_ * R).eval();
      }
    }
    else {
      if (current_step_%2 != 0) {
        ++current_step_;
        R = (R_yaw_ * R).eval();
      }
    }
    com.noalias() += contact_position[1];
    com.noalias() -= R * com_to_contact_position_local_[1];
    com.noalias() += contact_position[2];
    com.noalias() -= R * com_to_contact_position_local_[2];
    com.array() /= 2.0;
    contact_position[0].noalias() = com + R * (com_to_contact_position_local_[0] - 0.5 * step_length_);
    contact_position[3].noalias() = com + R * (com_to_contact_position_local_[3] - 0.5 * step_length_);
  }
  else {
    return false;
  } 
  com_ref_.clear();
  com_ref_.push_back(com);
  contact_position_ref_.clear();
  contact_position_ref_.push_back(contact_position);
  R_.clear();
  R_.push_back(R);
  if (enable_stance_phase_) {
    for (int step=current_step_; step<=planning_steps+current_step_; ++step) {
      if (step == 0) {
        // do nothing
      }
      else if (current_step_ == 0 && step == 1) {
        R = (R_yaw_ * R).eval();
        if (enable_raibert_heuristic_) {
          com.noalias() += 0.5 * R * step_length_;
        }
        else {
          com.noalias() += 0.25 * R * step_length_;
        }
        contact_position[1].noalias() = com + R * com_to_contact_position_local_[1];
        contact_position[2].noalias() = com + R * com_to_contact_position_local_[2];
      }
      else if (step%4 == 1) {
        R = (R_yaw_ * R).eval();
        com.noalias() += 0.5 * R * step_length_;
        contact_position[1].noalias() = com + R * com_to_contact_position_local_[1];
        contact_position[2].noalias() = com + R * com_to_contact_position_local_[2];
      }
      else if (step%4 == 3) {
        R = (R_yaw_ * R).eval();
        com.noalias() += 0.5 * R * step_length_;
        contact_position[0].noalias() = com + R * com_to_contact_position_local_[0];
        contact_position[3].noalias() = com + R * com_to_contact_position_local_[3];
      }
      com_ref_.push_back(com);
      contact_position_ref_.push_back(contact_position);
      R_.push_back(R);
    }
  }
  else {
    for (int step=current_step_; step<=planning_steps+current_step_; ++step) {
      if (step == 0) {
        // do nothing
      }
      else if (current_step_ == 0 && step == 1) {
        R = (R_yaw_ * R).eval();
        if (enable_raibert_heuristic_) {
          com.noalias() += 0.5 * R * step_length_;
        }
        else {
          com.noalias() += 0.25 * R * step_length_;
        }
        contact_position[1].noalias() = com + R * com_to_contact_position_local_[1];
        contact_position[2].noalias() = com + R * com_to_contact_position_local_[2];
      }
      else if (step%2 == 1) {
        R = (R_yaw_ * R).eval();
        com.noalias() += 0.5 * R * step_length_;
        contact_position[1].noalias() = com + R * com_to_contact_position_local_[1];
        contact_position[2].noalias() = com + R * com_to_contact_position_local_[2];
      }
      else {
        R = (R_yaw_ * R).eval();
        com.noalias() += 0.5 * R * step_length_;
        contact_position[0].noalias() = com + R * com_to_contact_position_local_[0];
        contact_position[3].noalias() = com + R * com_to_contact_position_local_[3];
      }
      com_ref_.push_back(com);
      contact_position_ref_.push_back(contact_position);
      R_.push_back(R);
    }
  }
  return true;
}


const aligned_vector<SE3>& TrotFootStepPlanner::contactPlacements(const int step) const {
  return contact_placement_ref_[step];
}


const aligned_vector<aligned_vector<SE3>>& TrotFootStepPlanner::contactPlacements() const {
  return contact_placement_ref_;
}


const std::vector<Eigen::Vector3d>& TrotFootStepPlanner::contactPositions(const int step) const {
  return contact_position_ref_[step];
}


const std::vector<std::vector<Eigen::Vector3d>>& TrotFootStepPlanner::contactPositions() const {
  return contact_position_ref_;
}


const Eigen::Vector3d& TrotFootStepPlanner::CoM(const int step) const {
  return com_ref_[step];
}
  

const std::vector<Eigen::Vector3d>& TrotFootStepPlanner::CoM() const {
  return com_ref_;
}


const Eigen::Matrix3d& TrotFootStepPlanner::R(const int step) const {
  return R_[step];
}
  

const std::vector<Eigen::Matrix3d>& TrotFootStepPlanner::R() const {
  return R_;
}


void TrotFootStepPlanner::disp(std::ostream& os) const {
  std::cout << "Trot foot step planner:" << std::endl;
  std::cout << "current_step:" << current_step_ << std::endl;
  const int planning_steps = contact_position_ref_.size();
  for (int i=0; i<planning_steps; ++i) {
    std::cout << "contact position[" << i << "]: ["  
              << contact_position_ref_[i][0].transpose() << "], [" 
              << contact_position_ref_[i][1].transpose() << "], [" 
              << contact_position_ref_[i][2].transpose() << "], [" 
              << contact_position_ref_[i][3].transpose() << "]" << std::endl;
    std::cout << "CoM position[" << i << "]: ["   << com_ref_[i].transpose() << "]" << std::endl;
    std::cout << "R[" << i << "]: ["   << R_[i] << "]" << std::endl;
  }
}


std::ostream& operator<<(std::ostream& os, 
                         const TrotFootStepPlanner& planner) {
  planner.disp(os);
  return os;
}


std::ostream& operator<<(std::ostream& os, 
                         const std::shared_ptr<TrotFootStepPlanner>& planner) {
  planner->disp(os);
  return os;
}

} // namespace robotoc 