#include "robotoc/mpc/walking_foot_step_planner.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

WalkingFootStepPlanner::WalkingFootStepPlanner(const Robot& biped_robot)
  : robot_(biped_robot),
    L_foot_id_(biped_robot.surfaceContactFrames()[0]),
    R_foot_id_(biped_robot.surfaceContactFrames()[1]),
    left_to_right_leg_distance_(0),
    contact_placement_ref_(),
    com_ref_(),
    com_to_contact_position_local_(),
    step_length_(Eigen::Vector3d::Zero()),
    com_(Eigen::Vector3d::Zero()),
    R_yaw_(Eigen::Matrix3d::Identity()) {
}


WalkingFootStepPlanner::WalkingFootStepPlanner() {
}


WalkingFootStepPlanner::~WalkingFootStepPlanner() {
}


void WalkingFootStepPlanner::setGaitPattern(const Eigen::Vector3d& step_length, 
                                            const double yaw_rate) {
  step_length_ = step_length;
  R_yaw_ << std::cos(yaw_rate), -std::sin(yaw_rate), 0, 
            std::sin(yaw_rate), std::cos(yaw_rate),  0,
            0, 0, 1;
  R_yaw_half_ << std::cos(0.5*yaw_rate), -std::sin(0.5*yaw_rate), 0, 
                 std::sin(0.5*yaw_rate), std::cos(0.5*yaw_rate),  0,
                 0, 0, 1;
}


void WalkingFootStepPlanner::init(const Eigen::VectorXd& q) {
  const Eigen::Matrix3d R 
      = Eigen::Quaterniond(q.coeff(6), q.coeff(3), q.coeff(4), q.coeff(5)).toRotationMatrix();
  robot_.updateFrameKinematics(q);
  com_to_contact_position_local_ = { R.transpose() * (robot_.framePosition(L_foot_id_)-robot_.CoM()), 
                                     R.transpose() * (robot_.framePosition(R_foot_id_)-robot_.CoM()) };
  left_to_right_leg_distance_ = com_to_contact_position_local_[0].coeff(1)
                                  - com_to_contact_position_local_[1].coeff(1);
}


bool WalkingFootStepPlanner::plan(const Eigen::VectorXd& q,
                                  const ContactStatus& contact_status,
                                  const int planning_steps) {
  assert(planning_steps >= 0);
  Eigen::Matrix3d R 
      = Eigen::Quaterniond(q.coeff(6), q.coeff(3), q.coeff(4), q.coeff(5)).toRotationMatrix();
  robot_.updateFrameKinematics(q);
  aligned_vector<SE3> contact_placement, contact_placement_local;
  for (const auto frame : robot_.surfaceContactFrames()) {
    contact_placement.push_back(robot_.framePlacement(frame));
    contact_placement_local.emplace_back(
        R.transpose() * robot_.frameRotation(frame),
        R.transpose() * (robot_.framePosition(frame) - q.template head<3>()));
  }
  int current_step = 0;
  if (contact_status.isContactActive(0) && contact_status.isContactActive(1)) {
    current_step = 0;
  }
  else if (contact_status.isContactActive(0)) {
    current_step = 1;
    // L
    contact_placement_local[1] = contact_placement_local[0];
    contact_placement_local[1].translation().coeffRef(1) -= left_to_right_leg_distance_;
    contact_placement_local[1].translation().noalias() -= 0.5 * step_length_;
    contact_placement[1].translation().noalias() = R * contact_placement_local[1].translation() + q.template head<3>();
    contact_placement[1].rotation().noalias() = R;
  }
  else if (contact_status.isContactActive(1)) {
    current_step = 2;
    // R
    contact_placement_local[0] = contact_placement_local[1];
    contact_placement_local[0].translation().coeffRef(1) += left_to_right_leg_distance_;
    contact_placement_local[0].translation().noalias() -= 0.5 * step_length_;
    contact_placement[0].translation().noalias() = R * contact_placement_local[0].translation() + q.template head<3>();
    contact_placement[0].rotation().noalias() = R;
  }
  else {
    return false;
  } 
  Eigen::Vector3d com = Eigen::Vector3d::Zero();
  for (int i=0; i<2; ++i) {
    com.noalias() += contact_placement[i].translation();
    com.noalias() -= R * com_to_contact_position_local_[i];
  }
  com.array() *= 0.5;
  contact_placement_ref_.clear();
  contact_placement_ref_.push_back(contact_placement);
  com_ref_.clear();
  com_ref_.push_back(com);
  for (int step=current_step; step<=planning_steps+current_step; ++step) {
    R = (R_yaw_ * R).eval();
    if (step == 0) {
      contact_placement_ref_.push_back(contact_placement);
      com_ref_.push_back(com);
      contact_placement[1].translation().noalias() += 0.5 * R * step_length_;
      contact_placement[1].rotation() = R;
      com.noalias() += 0.25 * R * step_length_;
    }
    else if (step%2 != 0) {
      contact_placement[1].translation().noalias() += R * step_length_;
      contact_placement[1].rotation() = R;
      com.noalias() += 0.5 * R * step_length_;
    }
    else {
      contact_placement[0].translation().noalias() += R * step_length_;
      contact_placement[0].rotation() = R;
      com.noalias() += 0.5 * R * step_length_;
    }
    contact_placement_ref_.push_back(contact_placement);
    com_ref_.push_back(com);
  }
  return true;
}


const aligned_vector<SE3>& WalkingFootStepPlanner::contactPlacement(
    const int step) const {
  return contact_placement_ref_[step];
}


const Eigen::Vector3d& WalkingFootStepPlanner::com(const int step) const {
  return com_ref_[step];
}
  

void WalkingFootStepPlanner::disp(std::ostream& os) const {
  std::cout << "Walking foot step planner:" << std::endl;
  const int planning_steps = contact_placement_ref_.size();
  for (int i=0; i<planning_steps; ++i) {
    std::cout << "contact placement[" << i << "]: ["  
              << contact_placement_ref_[i][0] << "], [" 
              << contact_placement_ref_[i][1] << "]" << std::endl;
    std::cout << "CoM position[" << i << "]: ["   << com_ref_[i].transpose() << "]" << std::endl;
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