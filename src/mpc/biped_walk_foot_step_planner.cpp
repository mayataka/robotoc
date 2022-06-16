#include "robotoc/mpc/biped_walk_foot_step_planner.hpp"
#include "robotoc/utils/rotation.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

BipedWalkFootStepPlanner::BipedWalkFootStepPlanner(const Robot& biped_robot)
  : ContactPlannerBase(),
    robot_(biped_robot),
    raibert_heuristic_(),
    vcom_moving_window_filter_(),
    enable_raibert_heuristic_(false),
    L_foot_id_(biped_robot.surfaceContactFrames()[0]),
    R_foot_id_(biped_robot.surfaceContactFrames()[1]),
    current_step_(0),
    left_to_right_leg_distance_(0),
    foot_height_to_com_height_(0),
    contact_placement_ref_(),
    com_ref_(),
    R_(),
    vcom_(Eigen::Vector3d::Zero()),
    vcom_cmd_(Eigen::Vector3d::Zero()),
    step_length_(Eigen::Vector3d::Zero()),
    R_yaw_(Eigen::Matrix3d::Identity()),
    R_current_(Eigen::Matrix3d::Identity()),
    yaw_rate_cmd_(0),
    enable_double_support_phase_(false) {
  try {
    if (biped_robot.maxNumSurfaceContacts() < 2) {
      throw std::out_of_range(
          "invalid argument: robot is not a bipedal robot!\n robot.maxNumSurfaceContacts() must be larger than 2!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


BipedWalkFootStepPlanner::BipedWalkFootStepPlanner() {
}


BipedWalkFootStepPlanner::~BipedWalkFootStepPlanner() {
}


void BipedWalkFootStepPlanner::setGaitPattern(const Eigen::Vector3d& step_length, 
                                              const double step_yaw, 
                                              const bool enable_double_support_phase) {
  step_length_ = step_length;
  R_yaw_<< std::cos(step_yaw), -std::sin(step_yaw), 0, 
           std::sin(step_yaw), std::cos(step_yaw),  0,
           0, 0, 1;
  enable_double_support_phase_ = enable_double_support_phase;
  enable_raibert_heuristic_ = false;
}


void BipedWalkFootStepPlanner::setGaitPattern(const Eigen::Vector3d& vcom_cmd, 
                                              const double yaw_rate_cmd, 
                                              const double swing_time, 
                                              const double double_support_time, 
                                              const double gain) {
  try {
    if (swing_time <= 0) {
      throw std::out_of_range("invalid value: swing_time must be positive!");
    }
    if (double_support_time < 0) {
      throw std::out_of_range("invalid value: double_support_time must be non-negative!");
    }
    if (gain <= 0.0) {
      throw std::out_of_range("invalid argument: gain must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  const double period = 2.0 * (swing_time + double_support_time);
  raibert_heuristic_.setParameters(period, gain);
  vcom_moving_window_filter_.setParameters(period, 0.1*period);
  vcom_cmd_ = vcom_cmd;
  const double step_yaw = yaw_rate_cmd * swing_time;
  R_yaw_<< std::cos(step_yaw), -std::sin(step_yaw), 0, 
           std::sin(step_yaw),  std::cos(step_yaw), 0,
           0, 0, 1;
  yaw_rate_cmd_ = yaw_rate_cmd;
  enable_double_support_phase_ = (double_support_time > 0.0);
  enable_raibert_heuristic_ = true;
}


void BipedWalkFootStepPlanner::init(const Eigen::VectorXd& q) {
  Eigen::Matrix3d R = rotation::toRotationMatrix(q.template segment<4>(3));
  rotation::projectRotationMatrix(R, rotation::ProjectionAxis::Z);
  robot_.updateFrameKinematics(q);
  std::vector<Eigen::Vector3d> contact_position_local 
      = { R.transpose() * (robot_.framePosition(L_foot_id_)-q.template head<3>()),
          R.transpose() * (robot_.framePosition(R_foot_id_)-q.template head<3>()) };
  left_to_right_leg_distance_ = contact_position_local[0].coeff(1)
                                  - contact_position_local[1].coeff(1);
  const double foot_height = 0.5 * (robot_.framePosition(L_foot_id_).coeff(2)
                                      + robot_.framePosition(R_foot_id_).coeff(2));
  foot_height_to_com_height_ = robot_.CoM().coeff(2) - foot_height;
  contact_position_ref_.clear();
  com_ref_.clear(),
  R_.clear();
  R_.push_back(R);
  vcom_moving_window_filter_.clear();
  current_step_ = 0;
}


bool BipedWalkFootStepPlanner::plan(const double t, const Eigen::VectorXd& q,
                                    const Eigen::VectorXd& v,
                                    const ContactStatus& contact_status,
                                    const int planning_steps) {
  assert(planning_steps >= 0);
  if (enable_raibert_heuristic_) {
    R_current_ = rotation::toRotationMatrix(q.template segment<4>(3));
    vcom_ = R_current_.transpose() * v.template head<3>();
    vcom_moving_window_filter_.push_back(t, vcom_.template head<2>());
    const Eigen::Vector2d& vcom_avg = vcom_moving_window_filter_.average();
    raibert_heuristic_.planStepLength(vcom_avg, vcom_cmd_.template head<2>(), 
                                      yaw_rate_cmd_);
    step_length_ = raibert_heuristic_.stepLength();
  }
  robot_.updateFrameKinematics(q);
  aligned_vector<SE3> contact_placement;
  for (const auto frame : robot_.surfaceContactFrames()) {
    contact_placement.push_back(robot_.framePlacement(frame));
  }
  Eigen::Vector3d com = Eigen::Vector3d::Zero();
  Eigen::Matrix3d R = R_.front();
  std::vector<Eigen::Vector3d> contact_position_local;
  for (const auto& e : contact_placement) {
    contact_position_local.push_back(
        R.transpose() * (e.translation()-q.template head<3>()));
  }
  if (contact_status.isContactActive(0) && contact_status.isContactActive(1)) {
    if (enable_double_support_phase_) {
      if (current_step_%2 != 0) {
        ++current_step_;
      }
    }
    else {
      current_step_ = 0;
    }
    com = 0.5 * (contact_placement[0].translation() + contact_placement[1].translation());
    com.coeffRef(2) += foot_height_to_com_height_;
  }
  else if (contact_status.isContactActive(0)) {
    if (enable_double_support_phase_) {
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
    contact_position_local[1].noalias() 
        = contact_position_local[0] - 0.5 * R_yaw_.transpose() * step_length_;
    contact_position_local[1].coeffRef(1) -= left_to_right_leg_distance_;
    contact_placement[1].translation().noalias()  
        = R * contact_position_local[1] + q.template head<3>();
    contact_placement[1].rotation() = R;
    com = 0.5 * (contact_placement[0].translation() + contact_placement[1].translation());
    com.coeffRef(2) += foot_height_to_com_height_;
  }
  else if (contact_status.isContactActive(1)) {
    if (enable_double_support_phase_) {
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
    contact_position_local[0].noalias() 
        = contact_position_local[1] - 0.5 * R_yaw_.transpose() * step_length_;
    contact_position_local[0].coeffRef(1) += left_to_right_leg_distance_;
    contact_placement[0].translation().noalias()  
        = R * contact_position_local[0] + q.template head<3>();
    contact_placement[0].rotation() = R;
    com = 0.5 * (contact_placement[0].translation() + contact_placement[1].translation());
    com.coeffRef(2) += foot_height_to_com_height_;
  }
  else {
    return false;
  } 
  com_ref_.clear();
  com_ref_.push_back(com);
  contact_placement_ref_.clear();
  contact_placement_ref_.push_back(contact_placement);
  R_.clear();
  R_.push_back(R);
  if (enable_double_support_phase_) {
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
        contact_placement[1].translation().noalias() += 0.5 * R * step_length_;
        contact_placement[1].rotation() = R;
      }
      else if (step%4 == 1) {
        R = (R_yaw_ * R).eval();
        com.noalias() += 0.5 * R * step_length_;
        contact_placement[1].translation().noalias() += R * step_length_;
        contact_placement[1].rotation() = R;
      }
      else if (step%4 == 3) {
        R = (R_yaw_ * R).eval();
        com.noalias() += 0.5 * R * step_length_;
        contact_placement[0].translation().noalias() += R * step_length_;
        contact_placement[0].rotation() = R;
      }
      contact_placement_ref_.push_back(contact_placement);
      com_ref_.push_back(com);
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
        contact_placement[1].translation().noalias() += 0.5 * R * step_length_;
        contact_placement[1].rotation() = R;
      }
      else if (step%2 == 1) {
        R = (R_yaw_ * R).eval();
        com.noalias() += 0.5 * R * step_length_;
        contact_placement[1].translation().noalias() += R * step_length_;
        contact_placement[1].rotation() = R;
      }
      else {
        R = (R_yaw_ * R).eval();
        com.noalias() += 0.5 * R * step_length_;
        contact_placement[0].translation().noalias() += R * step_length_;
        contact_placement[0].rotation() = R;
      }
      contact_placement_ref_.push_back(contact_placement);
      com_ref_.push_back(com);
      R_.push_back(R);
    }
  }
  contact_position_ref_.clear();
  for (const auto& e : contact_placement_ref_) {
    contact_position_ref_.push_back(
        std::vector<Eigen::Vector3d>({e[0].translation(), e[1].translation()}));
  }
  return true;
}


const aligned_vector<SE3>& BipedWalkFootStepPlanner::contactPlacement(const int step) const {
  return contact_placement_ref_[step];
}


const aligned_vector<aligned_vector<SE3>>& BipedWalkFootStepPlanner::contactPlacement() const {
  return contact_placement_ref_;
}


const std::vector<Eigen::Vector3d>& BipedWalkFootStepPlanner::contactPosition(const int step) const {
  return contact_position_ref_[step];
}


const std::vector<std::vector<Eigen::Vector3d>>& BipedWalkFootStepPlanner::contactPosition() const {
  return contact_position_ref_;
}


const Eigen::Vector3d& BipedWalkFootStepPlanner::com(const int step) const {
  return com_ref_[step];
}


const std::vector<Eigen::Vector3d>& BipedWalkFootStepPlanner::com() const {
  return com_ref_;
}


const Eigen::Matrix3d& BipedWalkFootStepPlanner::R(const int step) const {
  return R_[step];
}
  

const std::vector<Eigen::Matrix3d>& BipedWalkFootStepPlanner::R() const {
  return R_;
}


void BipedWalkFootStepPlanner::disp(std::ostream& os) const {
  std::cout << "BipedWalk foot step planner:" << std::endl;
  std::cout << "current_step:" << current_step_ << std::endl;
  const int planning_steps = contact_placement_ref_.size();
  for (int i=0; i<planning_steps; ++i) {
    std::cout << "contact placement[" << i << "]: ["  
              << contact_placement_ref_[i][0] << "], [" 
              << contact_placement_ref_[i][1] << "]" << std::endl;
    std::cout << "CoM position[" << i << "]: ["   << com_ref_[i].transpose() << "]" << std::endl;
    std::cout << "R[" << i << "]: ["   << R_[i] << "]" << std::endl;
  }
}


std::ostream& operator<<(std::ostream& os, 
                         const BipedWalkFootStepPlanner& planner) {
  planner.disp(os);
  return os;
}


std::ostream& operator<<(std::ostream& os, 
                         const std::shared_ptr<BipedWalkFootStepPlanner>& planner) {
  planner->disp(os);
  return os;
}

} // namespace robotoc 