#include "robotoc/mpc/jumping_foot_step_planner.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

JumpingFootStepPlanner::JumpingFootStepPlanner(const Robot& robot)
  : FootStepPlannerBase(),
    robot_(robot),
    contact_frames_(robot.contactFrames()),
    current_step_(0),
    contact_placement_ref_(),
    com_ref_(),
    com_to_contact_position_local_(),
    R_(),
    jump_length_(Eigen::Vector3d::Zero()),
    R_yaw_(Eigen::Matrix3d::Identity()) {
  try {
    if (robot.maxNumPointContacts() < 4 && robot.maxNumSurfaceContacts() < 2) {
      throw std::out_of_range(
          "invalid argument: robot is not a quadrupedal robot or a bipedal robot!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


JumpingFootStepPlanner::JumpingFootStepPlanner() {
}


JumpingFootStepPlanner::~JumpingFootStepPlanner() {
}


void JumpingFootStepPlanner::setJumpPattern(const Eigen::Vector3d& jump_length, 
                                            const double yaw_rate) {
  jump_length_ = jump_length;
  R_yaw_ << std::cos(yaw_rate), -std::sin(yaw_rate), 0, 
            std::sin(yaw_rate), std::cos(yaw_rate),  0,
            0, 0, 1;
}


void JumpingFootStepPlanner::init(const Eigen::VectorXd& q) {
  const Eigen::Matrix3d R 
      = Eigen::Quaterniond(q.coeff(6), q.coeff(3), q.coeff(4), q.coeff(5)).toRotationMatrix();
  robot_.updateFrameKinematics(q);
  aligned_vector<SE3> contact_placement;
  for (const auto frame : robot_.contactFrames()) {
    contact_placement.push_back(robot_.framePlacement(frame));
  }
  aligned_vector<SE3> contact_placement_local;
  for (const auto frame : robot_.contactFrames()) {
    contact_placement_local.emplace_back(
        R.transpose() * robot_.frameRotation(frame),
        R.transpose() * (robot_.framePosition(frame) - q.template head<3>()));
  }
  const Eigen::Matrix3d R_goal = R_yaw_ * R;
  const Eigen::Quaterniond quat_goal = Eigen::Quaterniond(R_goal);
  Eigen::VectorXd q_goal = q;
  q_goal.template head<3>().noalias() += R_yaw_ * jump_length_;
  q_goal.template segment<4>(3) = quat_goal.coeffs();
  robot_.updateFrameKinematics(q_goal);
  aligned_vector<SE3> contact_placement_goal;
  for (int i=0; i<robot_.contactFrames().size(); ++i) {
    contact_placement_goal.emplace_back(
        R_goal * contact_placement_local[i].rotation(),
        q.template head<3>() + R_goal * contact_placement_local[i].translation()
                             + R_yaw_ * jump_length_);
  }
  contact_placement_ref_.clear();
  contact_placement_ref_.push_back(contact_placement);
  contact_placement_ref_.push_back(contact_placement);
  contact_placement_ref_.push_back(contact_placement_goal);
  robot_.updateFrameKinematics(q);
  com_ref_.push_back(robot_.CoM());
  com_ref_.push_back(robot_.CoM());
  robot_.updateFrameKinematics(q_goal);
  com_ref_.push_back(robot_.CoM());
  R_.clear();
  R_.push_back(R);
  R_.push_back(R);
  R_.push_back(R_goal);
  current_step_ = 0;
}


bool JumpingFootStepPlanner::plan(const Eigen::VectorXd& q,
                                  const Eigen::VectorXd& v,
                                  const ContactStatus& contact_status,
                                  const int planning_steps) {
  assert(planning_steps >= 0);
  if (contact_status.hasActiveContacts()) {
    if (current_step_ == 1) {
      current_step_ = 2;
    }
    robot_.updateFrameKinematics(q);
    for (int i=0; i<contact_frames_.size(); ++i) {
      contact_placement_ref_[1][i] = robot_.framePlacement(contact_frames_[i]);
    }
    contact_placement_ref_[0] = contact_placement_ref_[1];
  }
  else {
    if (current_step_ == 0) {
      current_step_ = 1;
      contact_placement_ref_[1] = contact_placement_ref_[2];
      contact_placement_ref_.pop_back();
      com_ref_[1] = com_ref_[2];
      com_ref_.pop_back();
      R_[1] = R_[2];
      R_.pop_back();
    }
  }
  return true;
}


const aligned_vector<SE3>& JumpingFootStepPlanner::contactPlacement(const int step) const {
  return contact_placement_ref_[step];
}


const aligned_vector<aligned_vector<SE3>>& JumpingFootStepPlanner::contactPlacement() const {
  return contact_placement_ref_;
}


const std::vector<Eigen::Vector3d>& JumpingFootStepPlanner::contactPosition(const int step) const {
  return contact_position_ref_[step];
}


const std::vector<std::vector<Eigen::Vector3d>>& JumpingFootStepPlanner::contactPosition() const {
  return contact_position_ref_;
}


const Eigen::Vector3d& JumpingFootStepPlanner::com(const int step) const {
  return com_ref_[step];
}


const std::vector<Eigen::Vector3d>& JumpingFootStepPlanner::com() const {
  return com_ref_;
}


const Eigen::Matrix3d& JumpingFootStepPlanner::R(const int step) const {
  return R_[step];
}
  

const std::vector<Eigen::Matrix3d>& JumpingFootStepPlanner::R() const {
  return R_;
}


void JumpingFootStepPlanner::disp(std::ostream& os) const {
  std::cout << "Jumping foot step planner:" << std::endl;
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
                         const JumpingFootStepPlanner& planner) {
  planner.disp(os);
  return os;
}


std::ostream& operator<<(std::ostream& os, 
                         const std::shared_ptr<JumpingFootStepPlanner>& planner) {
  planner->disp(os);
  return os;
}

} // namespace robotoc 