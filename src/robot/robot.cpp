#include "robotoc/robot/robot.hpp"

#include <stdexcept>


namespace robotoc {

Robot::Robot(const std::string& path_to_urdf, 
             const BaseJointType& base_joint_type)
  : path_to_urdf_(path_to_urdf),
    model_(),
    impulse_model_(),
    data_(),
    impulse_data_(),
    contact_frame_names_(),
    contact_frames_(),
    contact_types_(),
    point_contacts_(),
    surface_contacts_(),
    fjoint_(),
    dimq_(0),
    dimv_(0),
    dimu_(0),
    dim_passive_(0),
    max_dimf_(0),
    max_num_contacts_(0),
    contact_inv_damping_(0),
    baumgarte_weights_({0, 0}),
    has_floating_base_(false),
    has_generalized_momentum_bias_(false),
    dimpulse_dv_(),
    generalized_momentum_bias_(),
    joint_effort_limit_(),
    joint_velocity_limit_(),
    lower_joint_position_limit_(),
    upper_joint_position_limit_() {
  switch (base_joint_type) {
    case BaseJointType::FloatingBase:
      pinocchio::urdf::buildModel(path_to_urdf, 
                                  pinocchio::JointModelFreeFlyer(), model_);
      dim_passive_ = 6;
      has_floating_base_ = true;
      break;
    default:
      pinocchio::urdf::buildModel(path_to_urdf, model_);
      dim_passive_ = 0;
      has_floating_base_ = false;
      break;
  }
  data_ = pinocchio::Data(model_);
  dimq_ = model_.nq;
  dimv_ = model_.nv;
  dimu_ = model_.nv - dim_passive_;
  generalized_momentum_bias_ = Eigen::VectorXd::Zero(dimv_);
  initializeJointLimits();
}


Robot::Robot(const std::string& path_to_urdf, 
             const BaseJointType& base_joint_type,
             const std::vector<int>& contact_frames,
             const std::vector<ContactType>& contact_types,
             const std::pair<double, double>& baumgarte_weights,
             const double contact_inv_damping)
  : Robot(path_to_urdf, base_joint_type) {
  if (baumgarte_weights.first < 0 || baumgarte_weights.second < 0) {
    throw std::out_of_range(
        "Invalid argument: baumgarte_weights must be non-negative!");
  }
  if (contact_frames.size() != contact_types.size()) {
    throw std::out_of_range(
        "Invalid argument: contact_frames.size() and contact_types.size() must be the same!");
  }
  if (contact_inv_damping < 0) {
    throw std::out_of_range(
        "Invalid argument: contact_inv_damping must be non-negative!");
  }
  baumgarte_weights_ = baumgarte_weights;
  impulse_model_ = model_;
  impulse_model_.gravity.linear().setZero();
  impulse_data_ = pinocchio::Data(impulse_model_);
  fjoint_ = pinocchio::container::aligned_vector<pinocchio::Force>(
                model_.joints.size(), pinocchio::Force::Zero());
  max_num_contacts_ = contact_frames.size();
  contact_frames_ = contact_frames;
  contact_frame_names_.clear();
  for (const auto e : contact_frames_) {
    contact_frame_names_.push_back(model_.frames[e].name);
  }
  contact_types_ = contact_types; 
  for (int i=0; i<contact_frames.size(); ++i) {
    switch (contact_types[i]) {
      case ContactType::PointContact:
        point_contacts_.push_back(PointContact(model_, contact_frames[i], 
                                               baumgarte_weights.first,
                                               baumgarte_weights.second));
        break;
      case ContactType::SurfaceContact:
        surface_contacts_.push_back(SurfaceContact(model_, contact_frames[i], 
                                                   baumgarte_weights.first,
                                                   baumgarte_weights.second));
        break;
      default:
        break;
    }
  }
  max_dimf_ = 3 * point_contacts_.size() + 6 * surface_contacts_.size();
  contact_inv_damping_ = contact_inv_damping;
  data_.JMinvJt.resize(max_dimf_, max_dimf_);
  data_.JMinvJt.setZero();
  data_.sDUiJt.resize(model_.nv, max_dimf_);
  data_.sDUiJt.setZero();
  dimpulse_dv_.resize(model_.nv, model_.nv);
  dimpulse_dv_.setZero();
}


Robot::Robot(const std::string& path_to_urdf, 
             const BaseJointType& base_joint_type,
             const std::vector<std::string>& contact_frame_names,
             const std::vector<ContactType>& contact_types,
             const std::pair<double, double>& baumgarte_weights,
             const double contact_inv_damping)
  : Robot(path_to_urdf, base_joint_type) {
  if (baumgarte_weights.first < 0 || baumgarte_weights.second < 0) {
    throw std::out_of_range(
        "Invalid argument: baumgarte_weights must be non-negative!");
  }
  if (contact_frame_names.size() != contact_types.size()) {
    throw std::out_of_range(
        "Invalid argument: contact_frame_names.size() and contact_types.size() must be the same!");
  }
  if (contact_inv_damping < 0) {
    throw std::out_of_range(
        "Invalid argument: contact_inv_damping must be non-negative!");
  }
  baumgarte_weights_ = baumgarte_weights;
  impulse_model_ = model_;
  impulse_model_.gravity.linear().setZero();
  impulse_data_ = pinocchio::Data(impulse_model_);
  fjoint_ = pinocchio::container::aligned_vector<pinocchio::Force>(
                model_.joints.size(), pinocchio::Force::Zero());
  max_num_contacts_ = contact_frame_names.size();
  contact_frames_.clear();
  for (const auto& e : contact_frame_names) {
    if (!model_.existFrame(e)) {
      throw std::invalid_argument(
          "Invalid argument: frame " + e + " does not exit!");
    }
    contact_frames_.push_back(model_.getFrameId(e));
  }
  contact_frame_names_ = contact_frame_names;
  contact_types_ = contact_types; 
  for (int i=0; i<contact_frames_.size(); ++i) {
    switch (contact_types[i]) {
      case ContactType::PointContact:
        point_contacts_.push_back(PointContact(model_, contact_frames_[i], 
                                               baumgarte_weights.first,
                                               baumgarte_weights.second));
        break;
      case ContactType::SurfaceContact:
        surface_contacts_.push_back(SurfaceContact(model_, contact_frames_[i], 
                                                   baumgarte_weights.first,
                                                   baumgarte_weights.second));
        break;
      default:
        break;
    }
  }
  max_dimf_ = 3 * point_contacts_.size() + 6 * surface_contacts_.size();
  contact_inv_damping_ = contact_inv_damping;
  data_.JMinvJt.resize(max_dimf_, max_dimf_);
  data_.JMinvJt.setZero();
  data_.sDUiJt.resize(model_.nv, max_dimf_);
  data_.sDUiJt.setZero();
  dimpulse_dv_.resize(model_.nv, model_.nv);
  dimpulse_dv_.setZero();
}


Robot::Robot(const std::string& path_to_urdf, 
             const BaseJointType& base_joint_type,
             const std::vector<int>& contact_frames,
             const std::vector<ContactType>& contact_types,
             const double time_step, const double contact_inv_damping)
  : Robot(path_to_urdf, base_joint_type, contact_frames, contact_types,
          std::make_pair(2.0/time_step, 1.0/(time_step*time_step)), 
          contact_inv_damping) {
}


Robot::Robot(const std::string& path_to_urdf, 
             const BaseJointType& base_joint_type,
             const std::vector<std::string>& contact_frame_names,
             const std::vector<ContactType>& contact_types,
             const double time_step, const double contact_inv_damping)
  : Robot(path_to_urdf, base_joint_type, contact_frame_names, contact_types,
          std::make_pair(2.0/time_step, 1.0/(time_step*time_step)), 
          contact_inv_damping) {
}


Robot::Robot()
  : path_to_urdf_(),
    model_(),
    impulse_model_(),
    data_(),
    impulse_data_(),
    contact_frames_(),
    contact_frame_names_(),
    contact_types_(),
    point_contacts_(),
    surface_contacts_(),
    fjoint_(),
    dimq_(0),
    dimv_(0),
    dimu_(0),
    dim_passive_(0),
    contact_inv_damping_(0),
    baumgarte_weights_({0, 0}),
    max_dimf_(0),
    max_num_contacts_(0),
    has_floating_base_(false),
    has_generalized_momentum_bias_(false),
    dimpulse_dv_(),
    generalized_momentum_bias_(),
    joint_effort_limit_(),
    joint_velocity_limit_(),
    lower_joint_position_limit_(),
    upper_joint_position_limit_() {
}


Robot::~Robot() {
}


Robot Robot::clone() const {
  BaseJointType base_joint_type;
  if (hasFloatingBase()) {
    base_joint_type = BaseJointType::FloatingBase;
  }
  else {
    base_joint_type = BaseJointType::FixedBase;
  }
  auto robot = Robot(path_to_urdf_, base_joint_type, 
                     contact_frame_names_, contact_types_, 
                     baumgarte_weights_, contact_inv_damping_);
  robot.setGeneralizedMomentumBias(generalizedMomentumBias());
  robot.setJointEffortLimit(jointEffortLimit());
  robot.setJointVelocityLimit(jointVelocityLimit());
  robot.setLowerJointPositionLimit(lowerJointPositionLimit());
  robot.setUpperJointPositionLimit(upperJointPositionLimit());
  return robot;
}


void Robot::initializeJointLimits() {
  const int dim_joint = model_.nv - dim_passive_;
  joint_effort_limit_.resize(dim_joint);
  joint_velocity_limit_.resize(dim_joint);
  lower_joint_position_limit_.resize(dim_joint);
  upper_joint_position_limit_.resize(dim_joint);
  joint_effort_limit_ = model_.effortLimit.tail(dim_joint);
  joint_velocity_limit_ = model_.velocityLimit.tail(dim_joint);
  lower_joint_position_limit_ = model_.lowerPositionLimit.tail(dim_joint);
  upper_joint_position_limit_ = model_.upperPositionLimit.tail(dim_joint);
}


void Robot::setJointEffortLimit(const Eigen::VectorXd& joint_effort_limit) {
  if (joint_effort_limit_.size() != joint_effort_limit.size()) {
    throw std::invalid_argument(
        "Invalid argument: joint_effort_limit.size() must be " 
        + std::to_string(joint_effort_limit_.size()));
  }
  joint_effort_limit_ = joint_effort_limit;
}


void Robot::setJointVelocityLimit(const Eigen::VectorXd& joint_velocity_limit) {
  if (joint_velocity_limit_.size() != joint_velocity_limit.size()) {
    throw std::invalid_argument(
        "Invalid argument: joint_velocity_limit.size() must be " 
        + std::to_string(joint_velocity_limit_.size()));
  }
  joint_velocity_limit_ = joint_velocity_limit;
}


void Robot::setLowerJointPositionLimit(
    const Eigen::VectorXd& lower_joint_position_limit) {
  if (lower_joint_position_limit_.size() != lower_joint_position_limit.size()) {
    throw std::invalid_argument(
        "Invalid argument: lower_joint_position_limit.size() must be " 
        + std::to_string(lower_joint_position_limit_.size()));
  }
  lower_joint_position_limit_ = lower_joint_position_limit;
}


void Robot::setUpperJointPositionLimit(
    const Eigen::VectorXd& upper_joint_position_limit) {
  if (upper_joint_position_limit_.size() != upper_joint_position_limit.size()) {
    throw std::invalid_argument(
        "Invalid argument: upper_joint_position_limit.size() must be " 
        + std::to_string(upper_joint_position_limit_.size()));
  }
  upper_joint_position_limit_ = upper_joint_position_limit;
}


void Robot::disp(std::ostream& os) const {
  os << "robot model:" << std::endl;
  os << "  name: " << model_.name << std::endl;
  if (has_floating_base_) {
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