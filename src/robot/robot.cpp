#include "robotoc/robot/robot.hpp"

#include <stdexcept>


namespace robotoc {

Robot::Robot(const RobotModelInfo& info)
  : info_(info),
    model_(),
    impulse_model_(),
    data_(),
    impulse_data_(),
    fjoint_(),
    dimpulse_dv_(),
    point_contacts_(),
    surface_contacts_(),
    dimq_(0),
    dimv_(0),
    dimu_(0),
    dim_passive_(0),
    max_dimf_(0),
    max_num_contacts_(0),
    properties_(),
    joint_effort_limit_(),
    joint_velocity_limit_(),
    lower_joint_position_limit_(),
    upper_joint_position_limit_() {
  switch (info.base_joint_type) {
    case BaseJointType::FloatingBase:
      pinocchio::urdf::buildModel(info.urdf_path, 
                                  pinocchio::JointModelFreeFlyer(), model_);
      dim_passive_ = 6;
      break;
    case BaseJointType::FixedBase:
      pinocchio::urdf::buildModel(info.urdf_path, model_);
      dim_passive_ = 0;
      break;
    default:
      std::runtime_error("[Robot] invalid argument: invalid base joint type");
      break;
  }
  impulse_model_ = model_;
  impulse_model_.gravity.linear().setZero();
  data_ = pinocchio::Data(model_);
  impulse_data_ = pinocchio::Data(impulse_model_);
  fjoint_ = pinocchio::container::aligned_vector<pinocchio::Force>(
                model_.joints.size(), pinocchio::Force::Zero());
  point_contacts_.clear();
  for (const auto& e : info.point_contacts) {
    point_contacts_.push_back(PointContact(model_, e));
  }
  surface_contacts_.clear();
  for (const auto& e : info.surface_contacts) {
    surface_contacts_.push_back(SurfaceContact(model_, e));
  }
  dimq_ = model_.nq;
  dimv_ = model_.nv;
  dimu_ = model_.nv - dim_passive_;
  max_dimf_ = 3 * point_contacts_.size() + 6 * surface_contacts_.size();
  max_num_contacts_ = point_contacts_.size() + surface_contacts_.size();
  data_.JMinvJt.resize(max_dimf_, max_dimf_);
  data_.JMinvJt.setZero();
  data_.sDUiJt.resize(model_.nv, max_dimf_);
  data_.sDUiJt.setZero();
  impulse_data_.JMinvJt.resize(max_dimf_, max_dimf_);
  impulse_data_.JMinvJt.setZero();
  impulse_data_.sDUiJt.resize(model_.nv, max_dimf_);
  impulse_data_.sDUiJt.setZero();
  dimpulse_dv_.resize(model_.nv, model_.nv);
  dimpulse_dv_.setZero();
  initializeJointLimits();
}


Robot::Robot()
  : info_(),
    model_(),
    impulse_model_(),
    data_(),
    impulse_data_(),
    fjoint_(),
    dimpulse_dv_(),
    point_contacts_(),
    surface_contacts_(),
    dimq_(0),
    dimv_(0),
    dimu_(0),
    dim_passive_(0),
    max_dimf_(0),
    max_num_contacts_(0),
    properties_(),
    joint_effort_limit_(),
    joint_velocity_limit_(),
    lower_joint_position_limit_(),
    upper_joint_position_limit_() {
}


void Robot::initializeJointLimits() {
  const int njoints = model_.nv - dim_passive_;
  joint_effort_limit_.resize(njoints);
  joint_velocity_limit_.resize(njoints);
  lower_joint_position_limit_.resize(njoints);
  upper_joint_position_limit_.resize(njoints);
  joint_effort_limit_ = model_.effortLimit.tail(njoints);
  joint_velocity_limit_ = model_.velocityLimit.tail(njoints);
  lower_joint_position_limit_ = model_.lowerPositionLimit.tail(njoints);
  upper_joint_position_limit_ = model_.upperPositionLimit.tail(njoints);
}


void Robot::setJointEffortLimit(const Eigen::VectorXd& joint_effort_limit) {
  if (joint_effort_limit_.size() != joint_effort_limit.size()) {
    throw std::invalid_argument(
        "[Robot] invalid argument: joint_effort_limit.size() must be " 
        + std::to_string(joint_effort_limit_.size()));
  }
  joint_effort_limit_ = joint_effort_limit;
}


void Robot::setJointVelocityLimit(const Eigen::VectorXd& joint_velocity_limit) {
  if (joint_velocity_limit_.size() != joint_velocity_limit.size()) {
    throw std::invalid_argument(
        "[Robot] invalid argument: joint_velocity_limit.size() must be " 
        + std::to_string(joint_velocity_limit_.size()));
  }
  joint_velocity_limit_ = joint_velocity_limit;
}


void Robot::setLowerJointPositionLimit(
    const Eigen::VectorXd& lower_joint_position_limit) {
  if (lower_joint_position_limit_.size() != lower_joint_position_limit.size()) {
    throw std::invalid_argument(
        "[Robot] invalid argument: lower_joint_position_limit.size() must be " 
        + std::to_string(lower_joint_position_limit_.size()));
  }
  lower_joint_position_limit_ = lower_joint_position_limit;
}


void Robot::setUpperJointPositionLimit(
    const Eigen::VectorXd& upper_joint_position_limit) {
  if (upper_joint_position_limit_.size() != upper_joint_position_limit.size()) {
    throw std::invalid_argument(
        "[Robot] invalid argument: upper_joint_position_limit.size() must be " 
        + std::to_string(upper_joint_position_limit_.size()));
  }
  upper_joint_position_limit_ = upper_joint_position_limit;
}


void Robot::disp(std::ostream& os) const {
  os << "Robot:" << std::endl;
  os << "  name: " << model_.name << std::endl;
  if (info_.base_joint_type == BaseJointType::FloatingBase) {
    os << "  base joint: floating base" << std::endl;
  }
  else {
    os << "  base joint: fixed base" << std::endl;
  }
  os << "  contacts:";
  for (const auto& e : point_contacts_) {
    os << e << std::endl;
  }
  for (const auto& e : surface_contacts_) {
    os << e << std::endl;
  }
  os << "  dimq = " << dimq_ << ", ";
  os << "  dimv = " << dimv_ << ", ";
  os << "  dimu = " << dimu_ << ", ";
  os << "  dim_passive = " << dim_passive_ << std::endl;
  os << std::endl;
  os << "  frames:" << std::endl;
  for (int i=0; i<model_.nframes; ++i) {
    os << "    frame " << i << std::endl;
    os << "      name: " << model_.frames[i].name << std::endl;
    os << "      parent joint id: " << model_.frames[i].parent << std::endl;
    os << std::endl;
  }
  os << "  joints:" << std::endl;
  for (int i=0; i<model_.njoints; ++i) {
    os << "    joint " << i << std::endl;
    os << "      name: " << model_.names[i] << std::endl;
    os << model_.joints[i] << std::endl;
  }
  os << "  effort limit = [" << joint_effort_limit_.transpose() << "]" 
            << std::endl;
  os << "  velocity limit = [" << joint_velocity_limit_.transpose() << "]"
            << std::endl;
  os << "  lowerl position limit = [" << lower_joint_position_limit_.transpose() 
            << "]" << std::endl;
  os << "  upper position limit = [" << upper_joint_position_limit_.transpose() 
            << "]" << std::flush; 
}


std::ostream& operator<<(std::ostream& os, const Robot& robot) {
  robot.disp(os);
  return os;
}

} // namespace robotoc 