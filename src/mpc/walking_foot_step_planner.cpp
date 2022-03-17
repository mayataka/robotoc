#include "robotoc/mpc/walking_foot_step_planner.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

WalkingFootStepPlanner::WalkingFootStepPlanner(const Robot& biped_robot)
  : FootStepPlannerBase(),
    robot_(biped_robot),
    L_foot_id_(biped_robot.surfaceContactFrames()[0]),
    R_foot_id_(biped_robot.surfaceContactFrames()[1]),
    current_step_(0),
    left_to_right_leg_distance_(0),
    foot_height_to_com_height_(0),
    contact_placement_ref_(),
    com_ref_(),
    R_(),
    step_length_(Eigen::Vector3d::Zero()),
    R_yaw_(Eigen::Matrix3d::Identity()),
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


WalkingFootStepPlanner::WalkingFootStepPlanner() {
}


WalkingFootStepPlanner::~WalkingFootStepPlanner() {
}


void WalkingFootStepPlanner::setGaitPattern(const Eigen::Vector3d& step_length, 
                                            const double yaw_rate,
                                            const bool enable_double_support_phase) {
  step_length_ = step_length;
  R_yaw_ << std::cos(yaw_rate), -std::sin(yaw_rate), 0, 
            std::sin(yaw_rate), std::cos(yaw_rate),  0,
            0, 0, 1;
  enable_double_support_phase_ = enable_double_support_phase;
}


void WalkingFootStepPlanner::init(const Eigen::VectorXd& q) {
  Eigen::Matrix3d R = Eigen::Quaterniond(q.coeff(6), q.coeff(3), q.coeff(4), q.coeff(5)).toRotationMatrix();
  R.coeffRef(0, 0) = 1.0;
  R.coeffRef(0, 1) = 0.0;
  R.coeffRef(0, 2) = 0.0;
  R.coeffRef(1, 0) = 0.0;
  R.coeffRef(2, 0) = 0.0;
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
  current_step_ = 0;
}


bool WalkingFootStepPlanner::plan(const Eigen::VectorXd& q,
                                  const Eigen::VectorXd& v,
                                  const ContactStatus& contact_status,
                                  const int planning_steps) {
  assert(planning_steps >= 0);
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
        com.noalias() += 0.25 * R * step_length_;
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
        com.noalias() += 0.25 * R * step_length_;
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


const aligned_vector<SE3>& WalkingFootStepPlanner::contactPlacement(const int step) const {
  return contact_placement_ref_[step];
}


const aligned_vector<aligned_vector<SE3>>& WalkingFootStepPlanner::contactPlacement() const {
  return contact_placement_ref_;
}


const std::vector<Eigen::Vector3d>& WalkingFootStepPlanner::contactPosition(const int step) const {
  return contact_position_ref_[step];
}


const std::vector<std::vector<Eigen::Vector3d>>& WalkingFootStepPlanner::contactPosition() const {
  return contact_position_ref_;
}


const Eigen::Vector3d& WalkingFootStepPlanner::com(const int step) const {
  return com_ref_[step];
}


const std::vector<Eigen::Vector3d>& WalkingFootStepPlanner::com() const {
  return com_ref_;
}


const Eigen::Matrix3d& WalkingFootStepPlanner::R(const int step) const {
  return R_[step];
}
  

const std::vector<Eigen::Matrix3d>& WalkingFootStepPlanner::R() const {
  return R_;
}


void WalkingFootStepPlanner::disp(std::ostream& os) const {
  std::cout << "Walking foot step planner:" << std::endl;
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
                         const WalkingFootStepPlanner& planner) {
  planner.disp(os);
  return os;
}


std::ostream& operator<<(std::ostream& os, 
                         const std::shared_ptr<WalkingFootStepPlanner>& planner) {
  planner->disp(os);
  return os;
}

} // namespace robotoc 