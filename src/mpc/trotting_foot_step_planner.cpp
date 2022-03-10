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
    contact_position_ref_(),
    com_ref_(),
    com_to_contact_position_local_(),
    step_length_(Eigen::Vector3d::Zero()),
    com_(Eigen::Vector3d::Zero()),
    R_yaw_(Eigen::Matrix3d::Identity()) {
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


void TrottingFootStepPlanner::setGaitPattern(
    const Eigen::Vector3d& step_length, const double yaw_rate) {
  step_length_ = step_length;
  R_yaw_<< std::cos(yaw_rate), -std::sin(yaw_rate), 0, 
           std::sin(yaw_rate), std::cos(yaw_rate),  0,
           0, 0, 1;
}


void TrottingFootStepPlanner::init(const Eigen::VectorXd& q) {
  const Eigen::Matrix3d R 
      = Eigen::Quaterniond(q.coeff(6), q.coeff(3), q.coeff(4), q.coeff(5)).toRotationMatrix();
  robot_.updateFrameKinematics(q);
  com_to_contact_position_local_ = { R.transpose() * (robot_.framePosition(LF_foot_id_)-robot_.CoM()), 
                                     R.transpose() * (robot_.framePosition(LH_foot_id_)-robot_.CoM()),
                                     R.transpose() * (robot_.framePosition(RF_foot_id_)-robot_.CoM()),
                                     R.transpose() * (robot_.framePosition(RH_foot_id_)-robot_.CoM()) };
}


bool TrottingFootStepPlanner::plan(const Eigen::VectorXd& q,
                                   const ContactStatus& contact_status,
                                   const int planning_steps) {
  assert(planning_steps >= 0);
  Eigen::Matrix3d R 
      = Eigen::Quaterniond(q.coeff(6), q.coeff(3), q.coeff(4), q.coeff(5)).toRotationMatrix();
  robot_.updateFrameKinematics(q);
  std::vector<Eigen::Vector3d> contact_position, contact_position_local;
  for (const auto frame : robot_.pointContactFrames()) {
    contact_position.push_back(robot_.framePosition(frame));
    contact_position_local.push_back(
        R.transpose() * (robot_.framePosition(frame) - q.template head<3>()));
  }
  int current_step = 0;
  if (contact_status.isContactActive(0) && contact_status.isContactActive(1) 
      && contact_status.isContactActive(2) && contact_status.isContactActive(3)) {
    current_step = 0;
  }
  else if (contact_status.isContactActive(0) && contact_status.isContactActive(3)) {
    current_step = 1;
    // LH
    contact_position_local[1] << contact_position_local[3].coeff(0), // x : same as RH 
                                 contact_position_local[0].coeff(1), // y : same as LF
                                 contact_position_local[3].coeff(2); // z : same as RH
    contact_position_local[1].noalias() -= 0.5 * step_length_;
    contact_position[1].noalias() = R * contact_position_local[1] + q.template head<3>();
    // RF
    contact_position_local[2] << contact_position_local[0].coeff(0), // x : same as LF 
                                 contact_position_local[3].coeff(1), // y : same as RH
                                 contact_position_local[0].coeff(2); // z : same as LF
    contact_position_local[2].noalias() -= 0.5 * step_length_;
    contact_position[2].noalias() = R * contact_position_local[2] + q.template head<3>();
  }
  else if (contact_status.isContactActive(1) && contact_status.isContactActive(2)) {
    current_step = 2;
    // LF
    contact_position_local[0] << contact_position_local[2].coeff(0), // x : same as RF
                                 contact_position_local[1].coeff(1), // y : same as LH
                                 contact_position_local[2].coeff(2); // z : same as RF
    contact_position_local[0].noalias() -= 0.5 * step_length_; 
    contact_position[0].noalias() = R * contact_position_local[0] + q.template head<3>();
    // RH
    contact_position_local[3] << contact_position_local[1].coeff(0), // x : same as LH 
                                 contact_position_local[2].coeff(1), // y : same as RF
                                 contact_position_local[1].coeff(2); // z : same as LH
    contact_position_local[3].noalias() -= 0.5 * step_length_; 
    contact_position[3].noalias() = R * contact_position_local[3] + q.template head<3>();
  }
  else {
    return false;
  } 
  Eigen::Vector3d com = Eigen::Vector3d::Zero();
  for (int i=0; i<4; ++i) {
    com.noalias() += contact_position[i];
    com.noalias() -= R * com_to_contact_position_local_[i];
  }
  com.array() *= 0.25;
  contact_position_ref_.clear();
  contact_position_ref_.push_back(contact_position);
  com_ref_.clear();
  com_ref_.push_back(com);
  for (int step=current_step; step<=planning_steps+current_step; ++step) {
    R = (R_yaw_ * R).eval();
    if (step == 0) {
      contact_position_ref_.push_back(contact_position);
      com_ref_.push_back(com);
      contact_position[1].noalias() += 0.5 * R * step_length_;
      contact_position[2].noalias() += 0.5 * R * step_length_;
      com.noalias() += 0.25 * R * step_length_;
    }
    else if (step%2 != 0) {
      contact_position[1].noalias() += R * step_length_;
      contact_position[2].noalias() += R * step_length_;
      com.noalias() += 0.5 * R * step_length_;
    }
    else {
      contact_position[0].noalias() += R * step_length_;
      contact_position[3].noalias() += R * step_length_;
      com.noalias() += 0.5 * R * step_length_;
    }
    contact_position_ref_.push_back(contact_position);
    com_ref_.push_back(com);
  }
  return true;
}


const aligned_vector<SE3>& TrottingFootStepPlanner::contactPlacement(
    const int step) const {
  return contact_placement_ref_[step];
}


const std::vector<Eigen::Vector3d>& TrottingFootStepPlanner::contactPosition(
    const int step) const {
  return contact_position_ref_[step];
}


const Eigen::Vector3d& TrottingFootStepPlanner::com(const int step) const {
  return com_ref_[step];
}
  

void TrottingFootStepPlanner::disp(std::ostream& os) const {
  std::cout << "Trotting foot step planner:" << std::endl;
  const int planning_steps = contact_position_ref_.size();
  for (int i=0; i<planning_steps; ++i) {
    std::cout << "contact position[" << i << "]: ["  
              << contact_position_ref_[i][0].transpose() << "], [" 
              << contact_position_ref_[i][1].transpose() << "], [" 
              << contact_position_ref_[i][2].transpose() << "], [" 
              << contact_position_ref_[i][3].transpose() << "]" << std::endl;
    std::cout << "CoM position[" << i << "]: ["   << com_ref_[i].transpose() << "]" << std::endl;
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