#include "robotoc/mpc/trotting_foot_step_planner.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

TrottingFootStepPlanner::TrottingFootStepPlanner(const Robot& quadruped_robot)
  : FootStepPlannerBase(),
    robot_(quadruped_robot),
    LF_foot_id_(quadruped_robot.pointContactFrames()[0]),
    LH_foot_id_(quadruped_robot.pointContactFrames()[1]),
    RF_foot_id_(quadruped_robot.pointContactFrames()[2]),
    RH_foot_id_(quadruped_robot.pointContactFrames()[3]),
    current_step_(0),
    contact_position_ref_(),
    com_ref_(),
    R_(),
    com_to_contact_position_local_(),
    step_length_(Eigen::Vector3d::Zero()),
    R_yaw_(Eigen::Matrix3d::Identity()),
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


TrottingFootStepPlanner::TrottingFootStepPlanner() {
}


TrottingFootStepPlanner::~TrottingFootStepPlanner() {
}


void TrottingFootStepPlanner::setGaitPattern(const Eigen::Vector3d& step_length, 
                                             const double yaw_rate, 
                                             const bool enable_stance_phase) {
  step_length_ = step_length;
  R_yaw_<< std::cos(yaw_rate), -std::sin(yaw_rate), 0, 
           std::sin(yaw_rate), std::cos(yaw_rate),  0,
           0, 0, 1;
  enable_stance_phase_ = enable_stance_phase;
}


void TrottingFootStepPlanner::init(const Eigen::VectorXd& q) {
  Eigen::Matrix3d R = Eigen::Quaterniond(q.coeff(6), q.coeff(3), q.coeff(4), q.coeff(5)).toRotationMatrix();
  R.coeffRef(0, 0) = 1.0;
  R.coeffRef(0, 1) = 0.0;
  R.coeffRef(0, 2) = 0.0;
  R.coeffRef(1, 0) = 0.0;
  R.coeffRef(2, 0) = 0.0;
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
  current_step_ = 0;
}


bool TrottingFootStepPlanner::plan(const Eigen::VectorXd& q,
                                   const Eigen::VectorXd& v,
                                   const ContactStatus& contact_status,
                                   const int planning_steps) {
  assert(planning_steps >= 0);
  robot_.updateFrameKinematics(q);
  std::vector<Eigen::Vector3d> contact_position;
  for (const auto frame : robot_.pointContactFrames()) {
    contact_position.push_back(robot_.framePosition(frame));
  }
  // Eigen::Vector3d com = com_ref_.front();
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
        // com.template head<2>() = robot_.CoM().template head<2>();
        R = (R_yaw_ * R).eval();
      }
    }
    else {
      if (current_step_%2 != 1) {
        ++current_step_;
        // com.template head<2>() = robot_.CoM().template head<2>();
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
        // com.template head<2>() = robot_.CoM().template head<2>();
        R = (R_yaw_ * R).eval();
      }
    }
    else {
      if (current_step_%2 != 0) {
        ++current_step_;
        // com.template head<2>() = robot_.CoM().template head<2>();
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
        com.noalias() += 0.25 * R * step_length_;
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
        com.noalias() += 0.25 * R * step_length_;
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


const aligned_vector<SE3>& TrottingFootStepPlanner::contactPlacement(const int step) const {
  return contact_placement_ref_[step];
}


const aligned_vector<aligned_vector<SE3>>& TrottingFootStepPlanner::contactPlacement() const {
  return contact_placement_ref_;
}


const std::vector<Eigen::Vector3d>& TrottingFootStepPlanner::contactPosition(const int step) const {
  return contact_position_ref_[step];
}


const std::vector<std::vector<Eigen::Vector3d>>& TrottingFootStepPlanner::contactPosition() const {
  return contact_position_ref_;
}


const Eigen::Vector3d& TrottingFootStepPlanner::com(const int step) const {
  return com_ref_[step];
}
  

const std::vector<Eigen::Vector3d>& TrottingFootStepPlanner::com() const {
  return com_ref_;
}


const Eigen::Matrix3d& TrottingFootStepPlanner::R(const int step) const {
  return R_[step];
}
  

const std::vector<Eigen::Matrix3d>& TrottingFootStepPlanner::R() const {
  return R_;
}


void TrottingFootStepPlanner::disp(std::ostream& os) const {
  std::cout << "Trotting foot step planner:" << std::endl;
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
                         const TrottingFootStepPlanner& planner) {
  planner.disp(os);
  return os;
}


std::ostream& operator<<(std::ostream& os, 
                         const std::shared_ptr<TrottingFootStepPlanner>& planner) {
  planner->disp(os);
  return os;
}

} // namespace robotoc 