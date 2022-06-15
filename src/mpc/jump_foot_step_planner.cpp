#include "robotoc/mpc/jump_foot_step_planner.hpp"
#include "robotoc/utils/rotation.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

JumpFootStepPlanner::JumpFootStepPlanner(const Robot& robot)
  : ContactPlannerBase(),
    robot_(robot),
    contact_frames_(robot.contactFrames()),
    current_step_(0),
    contact_placement_ref_(),
    com_ref_(),
    com_to_contact_position_local_(),
    R_(),
    jump_length_(Eigen::Vector3d::Zero()),
    R_yaw_(Eigen::Matrix3d::Identity()),
    is_biped_(false) {
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
  if (robot.maxNumPointContacts() < 4 && robot.maxNumSurfaceContacts() >= 2) {
    is_biped_ = true;
  }
}


JumpFootStepPlanner::JumpFootStepPlanner() {
}


JumpFootStepPlanner::~JumpFootStepPlanner() {
}


void JumpFootStepPlanner::setJumpPattern(const Eigen::Vector3d& jump_length, 
                                         const double yaw_step) {
  jump_length_ = jump_length;
  R_yaw_ << std::cos(yaw_step), -std::sin(yaw_step), 0, 
            std::sin(yaw_step), std::cos(yaw_step),  0,
            0, 0, 1;
}


void JumpFootStepPlanner::init(const Eigen::VectorXd& q) {
  Eigen::Matrix3d R = rotation::toRotationMatrix(q.template segment<4>(3));
  rotation::projectRotationMatrix(R, rotation::ProjectionAxis::Z);
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
  if (!is_biped_) {
    for (auto& e : contact_placement_ref_) {
      for (auto& ee : e) {
        ee.rotation().setIdentity();
      }
    }
  }
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


bool JumpFootStepPlanner::plan(const double t, const Eigen::VectorXd& q,
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
      if (is_biped_) {
        contact_placement_ref_[1][i] 
            = robot_.framePlacement(contact_frames_[i]);
      }
      else {
        contact_placement_ref_[1][i].translation() 
            = robot_.framePosition(contact_frames_[i]);
      }
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


const aligned_vector<SE3>& JumpFootStepPlanner::contactPlacement(const int step) const {
  return contact_placement_ref_[step];
}


const aligned_vector<aligned_vector<SE3>>& JumpFootStepPlanner::contactPlacement() const {
  return contact_placement_ref_;
}


const std::vector<Eigen::Vector3d>& JumpFootStepPlanner::contactPosition(const int step) const {
  return contact_position_ref_[step];
}


const std::vector<std::vector<Eigen::Vector3d>>& JumpFootStepPlanner::contactPosition() const {
  return contact_position_ref_;
}


const Eigen::Vector3d& JumpFootStepPlanner::com(const int step) const {
  return com_ref_[step];
}


const std::vector<Eigen::Vector3d>& JumpFootStepPlanner::com() const {
  return com_ref_;
}


const Eigen::Matrix3d& JumpFootStepPlanner::R(const int step) const {
  return R_[step];
}
  

const std::vector<Eigen::Matrix3d>& JumpFootStepPlanner::R() const {
  return R_;
}


void JumpFootStepPlanner::disp(std::ostream& os) const {
  std::cout << "Jump foot step planner:" << std::endl;
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
                         const JumpFootStepPlanner& planner) {
  planner.disp(os);
  return os;
}


std::ostream& operator<<(std::ostream& os, 
                         const std::shared_ptr<JumpFootStepPlanner>& planner) {
  planner->disp(os);
  return os;
}

} // namespace robotoc 