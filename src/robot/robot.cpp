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


const Eigen::Vector3d& Robot::framePosition(const int frame_id) const {
  return data_.oMf[frame_id].translation();
}


const Eigen::Vector3d& Robot::framePosition(
    const std::string& frame_name) const {
  return framePosition(frameId(frame_name));
}


const Eigen::Matrix3d& Robot::frameRotation(const int frame_id) const {
  return data_.oMf[frame_id].rotation();
}


const Eigen::Matrix3d& Robot::frameRotation(
    const std::string& frame_name) const {
  return frameRotation(frameId(frame_name));
}


const SE3& Robot::framePlacement(const int frame_id) const {
  return data_.oMf[frame_id];
}


const SE3& Robot::framePlacement(const std::string& frame_name) const {
  return framePlacement(frameId(frame_name));
}


const Eigen::Vector3d& Robot::CoM() const {
  return data_.com[0];
}


Eigen::Vector3d Robot::frameLinearVelocity(
    const int frame_id, const pinocchio::ReferenceFrame reference_frame) const {
  return pinocchio::getFrameVelocity(model_, data_, frame_id, reference_frame).linear();
}


Eigen::Vector3d Robot::frameLinearVelocity(
    const std::string& frame_name, 
    const pinocchio::ReferenceFrame reference_frame) const {
  return frameLinearVelocity(frameId(frame_name), reference_frame);
}


Eigen::Vector3d Robot::frameAngularVelocity(
    const int frame_id, const pinocchio::ReferenceFrame reference_frame) const {
  return pinocchio::getFrameVelocity(model_, data_, frame_id, reference_frame).angular();
}


Eigen::Vector3d Robot::frameAngularVelocity(
    const std::string& frame_name, 
    const pinocchio::ReferenceFrame reference_frame) const {
  return frameAngularVelocity(frameId(frame_name), reference_frame);
}


Robot::Vector6d Robot::frameSpatialVelocity(
    const int frame_id, const pinocchio::ReferenceFrame reference_frame) const {
  return pinocchio::getFrameVelocity(model_, data_, frame_id, reference_frame).toVector();
}


Robot::Vector6d Robot::frameSpatialVelocity(
    const std::string& frame_name, 
    const pinocchio::ReferenceFrame reference_frame) const {
  return frameSpatialVelocity(frameId(frame_name), reference_frame);
}


const Eigen::Vector3d& Robot::CoMVelocity() const {
  return data_.vcom[0];
}


void Robot::setContactForces(const ContactStatus& contact_status, 
                             const std::vector<Vector6d>& f) {
  assert(f.size() == max_num_contacts_);
  const int num_point_contacts = point_contacts_.size();
  const int num_surface_contacts = surface_contacts_.size();
  for (int i=0; i<num_point_contacts; ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_[i].computeJointForceFromContactForce(
          f[i].template head<3>(), fjoint_);
    }
    else {
      point_contacts_[i].computeJointForceFromContactForce(
          Eigen::Vector3d::Zero(), fjoint_);
    }
  }
  for (int i=0; i<num_surface_contacts; ++i) {
    if (contact_status.isContactActive(i+num_point_contacts)) {
      surface_contacts_[i].computeJointForceFromContactWrench(
          f[i+num_point_contacts], fjoint_);
    }
    else {
      surface_contacts_[i].computeJointForceFromContactWrench(
          Vector6d::Zero(), fjoint_);
    }
  }
}


Eigen::VectorXd Robot::generateFeasibleConfiguration() const {
  Eigen::VectorXd q_min(dimq_), q_max(dimq_);
  if (info_.base_joint_type == BaseJointType::FloatingBase) {
    q_min.template head<7>() = - Eigen::VectorXd::Ones(7);
    q_max.template head<7>() = Eigen::VectorXd::Ones(7);
  }
  q_min.tail(dimu_) = lower_joint_position_limit_;
  q_max.tail(dimu_) = upper_joint_position_limit_;
  return pinocchio::randomConfiguration(model_, q_min, q_max);
}


int Robot::frameId(const std::string& frame_name) const {
  if (!model_.existFrame(frame_name)) {
    throw std::invalid_argument(
        "[Robot] invalid argument: frame '" + frame_name + "' does not exit!");
  }
  return model_.getFrameId(frame_name);
}


std::string Robot::frameName(const int frame_id) const {
  return  model_.frames[frame_id].name;
}


double Robot::totalMass() const {
  return pinocchio::computeTotalMass(model_);
}


double Robot::totalWeight() const {
  return (- pinocchio::computeTotalMass(model_) * model_.gravity981.coeff(2));
}


int Robot::dimq() const {
  return dimq_;
}


int Robot::dimv() const {
  return dimv_;
}


int Robot::dimu() const {
  return dimu_;
}


int Robot::max_dimf() const {
  return max_dimf_;
}


int Robot::dim_passive() const {
  return dim_passive_;
}


bool Robot::hasFloatingBase() const {
  return (info_.base_joint_type == BaseJointType::FloatingBase);
}


int Robot::maxNumContacts() const {
  return (maxNumPointContacts() + maxNumSurfaceContacts());
}


int Robot::maxNumPointContacts() const {
  return point_contacts_.size();
}


int Robot::maxNumSurfaceContacts() const {
  return surface_contacts_.size();
}


ContactType Robot::contactType(const int contact_index) const {
  assert(contact_index >= 0);
  assert(contact_index < max_num_contacts_);
  if (contact_index < maxNumPointContacts()) {
    return ContactType::PointContact;
  }
  else {
    return ContactType::SurfaceContact;
  }
}


std::vector<ContactType> Robot::contactTypes() const {
  std::vector<ContactType> contact_types;
  for (const auto& e : point_contacts_) {
    contact_types.push_back(ContactType::PointContact);
  }
  for (const auto& e : surface_contacts_) {
    contact_types.push_back(ContactType::SurfaceContact);
  }
  return contact_types;
}


std::vector<int> Robot::contactFrames() const {
  std::vector<int> contact_frames;
  for (const auto& e : point_contacts_) {
    contact_frames.push_back(e.contactFrameId());
  }
  for (const auto& e : surface_contacts_) {
    contact_frames.push_back(e.contactFrameId());
  }
  return contact_frames;
}


std::vector<std::string> Robot::contactFrameNames() const {
  std::vector<std::string> contact_frames;
  for (const auto& e : point_contacts_) {
    contact_frames.push_back(e.contactModelInfo().frame);
  }
  for (const auto& e : surface_contacts_) {
    contact_frames.push_back(e.contactModelInfo().frame);
  }
  return contact_frames;
}


std::vector<int> Robot::pointContactFrames() const {
  std::vector<int> contact_frames;
  for (const auto& e : point_contacts_) {
    contact_frames.push_back(e.contactFrameId());
  }
  return contact_frames;
}


std::vector<std::string> Robot::pointContactFrameNames() const {
  std::vector<std::string> contact_frames;
  for (const auto& e : point_contacts_) {
    contact_frames.push_back(e.contactModelInfo().frame);
  }
  return contact_frames;
}


std::vector<int> Robot::surfaceContactFrames() const {
  std::vector<int> contact_frames;
  for (const auto& e : surface_contacts_) {
    contact_frames.push_back(e.contactFrameId());
  }
  return contact_frames;
}


std::vector<std::string> Robot::surfaceContactFrameNames() const {
  std::vector<std::string> contact_frames;
  for (const auto& e : surface_contacts_) {
    contact_frames.push_back(e.contactModelInfo().frame);
  }
  return contact_frames;
}


ContactStatus Robot::createContactStatus() const {
  return ContactStatus(contactTypes(), contactFrameNames());
}


ImpulseStatus Robot::createImpulseStatus() const {
  return ImpulseStatus(contactTypes(), contactFrameNames());
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


const RobotModelInfo& Robot::robotModelInfo() const {
  return info_;
}


const RobotProperties& Robot::robotProperties() const {
  return properties_;
}


void Robot::setRobotProperties(const RobotProperties& properties) {
  properties_ = properties;
  properties_.has_generalized_momentum_bias = false;
  if (properties_.generalized_momentum_bias.size() == dimv_) {
    properties_.has_generalized_momentum_bias 
        = !(properties_.generalized_momentum_bias.isZero());
  } 
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