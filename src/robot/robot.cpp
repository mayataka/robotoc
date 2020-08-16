#include "idocp/robot/robot.hpp"

#include <limits>
#include <assert.h>


namespace idocp {

Robot::Robot(const std::string& urdf_file_name)
  : model_(),
    data_(model_),
    urdf_file_name_(urdf_file_name),
    point_contacts_(),
    floating_base_(),
    fjoint_(),
    dimq_(0),
    dimv_(0),
    dimJ_(0),
    max_dimf_(0),
    dimf_(0),
    num_active_contacts_(0),
    is_each_contact_active_(),
    joint_effort_limit_(),
    joint_velocity_limit_(),
    lower_joint_position_limit_(),
    upper_joint_position_limit_() {
  pinocchio::urdf::buildModel(urdf_file_name, model_);
  data_ = pinocchio::Data(model_);
  fjoint_ = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  floating_base_ = FloatingBase(model_);
  dimq_ = model_.nq;
  dimv_ = model_.nv;
  dimJ_ = model_.joints.size();
  initializeJointLimits();
}


Robot::Robot(const std::string& urdf_file_name, 
             const std::vector<int>& contact_frames, 
             const double baumgarte_weight_on_velocity, 
             const double baumgarte_weight_on_position) 
  : model_(),
    data_(model_),
    urdf_file_name_(urdf_file_name),
    point_contacts_(),
    floating_base_(),
    fjoint_(),
    dimq_(0),
    dimv_(0),
    dimJ_(0),
    max_dimf_(0),
    dimf_(0),
    num_active_contacts_(0),
    is_each_contact_active_(),
    joint_effort_limit_(),
    joint_velocity_limit_(),
    lower_joint_position_limit_(),
    upper_joint_position_limit_() {
  assert(baumgarte_weight_on_velocity >= 0);
  assert(baumgarte_weight_on_position >= 0);
  pinocchio::urdf::buildModel(urdf_file_name, model_);
  data_ = pinocchio::Data(model_);
  for (const auto& frame : contact_frames) {
    point_contacts_.push_back(PointContact(model_, frame, 
                                           baumgarte_weight_on_velocity,
                                           baumgarte_weight_on_position));
    is_each_contact_active_.push_back(false);
  }
  max_dimf_ = 3 * point_contacts_.size();
  fjoint_ = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  floating_base_ = FloatingBase(model_);
  dimq_ = model_.nq;
  dimv_ = model_.nv;
  dimJ_ = model_.joints.size();
  initializeJointLimits();
}


Robot::Robot()
  : model_(),
    data_(model_),
    urdf_file_name_(),
    point_contacts_(),
    floating_base_(),
    fjoint_(),
    dimq_(0),
    dimv_(0),
    dimJ_(0),
    max_dimf_(0),
    dimf_(0),
    num_active_contacts_(0),
    is_each_contact_active_(),
    joint_effort_limit_(),
    joint_velocity_limit_(),
    lower_joint_position_limit_(),
    upper_joint_position_limit_() {
}


Robot::~Robot() {
}


void Robot::buildRobotModelFromXML(const std::string& xml) {
  pinocchio::urdf::buildModelFromXML(xml, model_);
  data_ = pinocchio::Data(model_);
  fjoint_ = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  floating_base_ = FloatingBase(model_);
  dimq_ = model_.nq;
  dimv_ = model_.nv;
  dimJ_ = model_.joints.size();
  initializeJointLimits();
}


void Robot::buildRobotModelFromXML(const std::string& xml,
                                   const std::vector<int>& contact_frames, 
                                   const double baumgarte_weight_on_velocity, 
                                   const double baumgarte_weight_on_position) {
  pinocchio::urdf::buildModelFromXML(xml, model_);
  data_ = pinocchio::Data(model_);
  for (const auto& frame : contact_frames) {
    point_contacts_.push_back(PointContact(model_, frame, 
                                           baumgarte_weight_on_velocity,
                                           baumgarte_weight_on_position));
    is_each_contact_active_.push_back(false);
  }
  max_dimf_ = 3 * point_contacts_.size();
  fjoint_ = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  floating_base_ = FloatingBase(model_);
  dimq_ = model_.nq;
  dimv_ = model_.nv;
  dimJ_ = model_.joints.size();
  initializeJointLimits();
}


void Robot::setContactPoints(
    const std::vector<Eigen::Vector3d>& contact_points) {
  for (int i=0; i<point_contacts_.size(); ++i) {
    point_contacts_[i].resetContactPoint(contact_points[i]);
  }
}


void Robot::setContactPointsByCurrentKinematics() {
  for (int i=0; i<point_contacts_.size(); ++i) {
    point_contacts_[i].resetContactPointByCurrentKinematics(data_);
  }
}


void Robot::setContactStatus(const std::vector<bool>& is_each_contact_active) {
  assert(is_each_contact_active.size() == is_each_contact_active_.size());
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    is_each_contact_active_[i] = is_each_contact_active[i];
    if (is_each_contact_active[i]) {
      point_contacts_[i].activate();
      ++num_active_contacts;
    }
    else {
      point_contacts_[i].deactivate();
    }
  }
  num_active_contacts_ = num_active_contacts;
  dimf_ = 3 * num_active_contacts;
}


Eigen::VectorXd Robot::jointEffortLimit() const {
  return joint_effort_limit_;
}


Eigen::VectorXd Robot::jointVelocityLimit() const {
  return joint_velocity_limit_;
}

Eigen::VectorXd Robot::lowerJointPositionLimit() const {
  return lower_joint_position_limit_;
}


Eigen::VectorXd Robot::upperJointPositionLimit() const {
  return upper_joint_position_limit_;
}


void Robot::setJointEffortLimit(const Eigen::VectorXd& joint_effort_limit) {
  assert(joint_effort_limit_.size() == joint_effort_limit.size());
  joint_effort_limit_ = joint_effort_limit;
}


void Robot::setJointVelocityLimit(const Eigen::VectorXd& joint_velocity_limit) {
  assert(joint_velocity_limit_.size() == joint_velocity_limit.size());
  joint_velocity_limit_ = joint_velocity_limit;
}


void Robot::setLowerJointPositionLimit(
    const Eigen::VectorXd& lower_joint_position_limit) {
  assert(
      lower_joint_position_limit_.size() == lower_joint_position_limit.size());
  lower_joint_position_limit_ = lower_joint_position_limit;
}


void Robot::setUpperJointPositionLimit(
    const Eigen::VectorXd& upper_joint_position_limit) {
  assert(
      upper_joint_position_limit_.size() == upper_joint_position_limit.size());
  upper_joint_position_limit_ = upper_joint_position_limit;
}


int Robot::dimq() const {
  return dimq_;
}


int Robot::dimv() const {
  return dimv_;
}


int Robot::dimJ() const {
  return dimJ_;
}


int Robot::max_dimf() const {
  return max_dimf_;
}


int Robot::dimf() const {
  return dimf_;
}


int Robot::dim_passive() const {
  return floating_base_.dim_passive();
}


bool Robot::has_floating_base() const {
  return floating_base_.has_floating_base();
}


int Robot::max_point_contacts() const {
  return point_contacts_.size();
}


int Robot::num_active_point_contacts() const {
  return num_active_contacts_;
}


bool Robot::is_contact_active(const int contact_index) const {
  assert(contact_index >= 0);
  assert(contact_index < point_contacts_.size());
  return point_contacts_[contact_index].isActive();
}


void Robot::printRobotModel() const {
  for (int i=0; i<model_.njoints; ++i) {
    std::cout << "Info of joint " << i << std::endl;
    std::cout << "name: " << model_.names[i] << std::endl;
    std::cout << model_.joints[i] << std::endl;
  }
  std::cout << "effortLimit = [" << model_.effortLimit.transpose() << "]" 
            << std::endl;
  std::cout << "velocityLimit = [" << model_.velocityLimit.transpose() << "]"
            << std::endl;
  std::cout << "lowerPositionLimit = [" << model_.lowerPositionLimit.transpose() 
            << "]" << std::endl;
  std::cout << "upperPositionLimit = [" << model_.upperPositionLimit.transpose() 
            << "]" << std::endl;
}


void Robot::initializeJointLimits() {
  const int dim_joint = model_.nv - floating_base_.dim_passive();
  joint_effort_limit_.resize(dim_joint);
  joint_velocity_limit_.resize(dim_joint);
  lower_joint_position_limit_.resize(dim_joint);
  upper_joint_position_limit_.resize(dim_joint);
  joint_effort_limit_ = model_.effortLimit.tail(dim_joint);
  joint_velocity_limit_ = model_.velocityLimit.tail(dim_joint);
  lower_joint_position_limit_ = model_.lowerPositionLimit.tail(dim_joint);
  upper_joint_position_limit_ = model_.upperPositionLimit.tail(dim_joint);
}

} // namespace idocp 