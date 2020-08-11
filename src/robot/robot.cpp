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


void Robot::integrateConfiguration(const Eigen::Ref<const Eigen::VectorXd>& v, 
                                   const double integration_length, 
                                   Eigen::Ref<Eigen::VectorXd> q) const {
  assert(v.size() == dimv_);
  assert(integration_length >= 0);
  assert(q.size() == dimq_);
  if (floating_base_.has_floating_base()) {
    q = pinocchio::integrate(model_, q, integration_length*v);
  }
  else {
    q.noalias() += integration_length * v;
  }
}


void Robot::subtractConfiguration(
    const Eigen::Ref<const Eigen::VectorXd>& q_plus, 
    const Eigen::Ref<const Eigen::VectorXd>& q_minus, 
    Eigen::Ref<Eigen::VectorXd> difference) const {
  assert(q_plus.size() == dimq_);
  assert(q_minus.size() == dimq_);
  assert(difference.size() == dimv_);
  if (floating_base_.has_floating_base()) {
    difference = pinocchio::difference(model_, q_minus, q_plus);
  }
  else {
    difference = q_plus - q_minus;
  }
}


void Robot::dIntegrateConfiguration(
    const Eigen::Ref<const Eigen::VectorXd>& q, 
    const Eigen::Ref<const Eigen::VectorXd>& v, const double integration_length,
    Eigen::Ref<Eigen::MatrixXd> dIntegrate_dq, 
    Eigen::Ref<Eigen::MatrixXd> dIntegrate_dv) const {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(integration_length >= 0);
  assert(dIntegrate_dq.rows() == dimv_);
  assert(dIntegrate_dq.cols() == dimv_);
  assert(dIntegrate_dv.rows() == dimv_);
  assert(dIntegrate_dv.cols() == dimv_);
  pinocchio::dIntegrate(model_, q, v, dIntegrate_dq, pinocchio::ARG0);
  pinocchio::dIntegrate(model_, q, v, dIntegrate_dv, pinocchio::ARG1);
}


void Robot::dSubtractdConfigurationPlus(
      const Eigen::Ref<const Eigen::VectorXd>& q_plus,
      const Eigen::Ref<const Eigen::VectorXd>& q_minus,
      Eigen::Ref<Eigen::MatrixXd> dSubtract_dqplus) const {
  assert(q_plus.size() == dimq_);
  assert(q_minus.size() == dimq_);
  assert(dSubtract_dqplus.rows() == dimv_);
  assert(dSubtract_dqplus.cols() == dimv_);
  pinocchio::dDifference(model_, q_minus, q_plus, dSubtract_dqplus, 
                         pinocchio::ARG1);
}


void Robot::dSubtractdConfigurationMinus(
      const Eigen::Ref<const Eigen::VectorXd>& q_plus,
      const Eigen::Ref<const Eigen::VectorXd>& q_minus,
      Eigen::Ref<Eigen::MatrixXd> dSubtract_dqminus) const {
  assert(q_plus.size() == dimq_);
  assert(q_minus.size() == dimq_);
  assert(dSubtract_dqminus.rows() == dimv_);
  assert(dSubtract_dqminus.cols() == dimv_);
  pinocchio::dDifference(model_, q_minus, q_plus, dSubtract_dqminus, 
                         pinocchio::ARG0);
}


void Robot::computeConfigurationJacobian(
    const Eigen::Ref<const Eigen::VectorXd>& q, 
    Eigen::Ref<Eigen::MatrixXd> J) const {
  assert(q.size() == dimq_);
  assert(J.rows() == dimq_);
  assert(J.cols() == dimv_);
  pinocchio::integrateCoeffWiseJacobian(model_, q, J);
}


void Robot::updateKinematics(const Eigen::Ref<const Eigen::VectorXd>& q, 
                             const Eigen::Ref<const Eigen::VectorXd>& v, 
                             const Eigen::Ref<const Eigen::VectorXd>& a) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  pinocchio::forwardKinematics(model_, data_, q, v, a);
  pinocchio::updateFramePlacements(model_, data_);
  pinocchio::computeForwardKinematicsDerivatives(model_, data_, q, v, a);
}


void Robot::computeBaumgarteResidual(
    Eigen::Ref<Eigen::VectorXd> baumgarte_residual) const {
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].computeBaumgarteResidual(
          model_, data_, baumgarte_residual.segment<3>(3*num_active_contacts));
      ++num_active_contacts;
    }
  }
}


void Robot::computeBaumgarteResidual(
    const double coeff, Eigen::Ref<Eigen::VectorXd> baumgarte_residual) const {
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].computeBaumgarteResidual(
          model_, data_, coeff, 
          baumgarte_residual.segment<3>(3*num_active_contacts));
      ++num_active_contacts;
    }
  }
}


void Robot::computeBaumgarteDerivatives(
    Eigen::Ref<Eigen::MatrixXd> baumgarte_partial_dq, 
    Eigen::Ref<Eigen::MatrixXd> baumgarte_partial_dv, 
    Eigen::Ref<Eigen::MatrixXd> baumgarte_partial_da) {
  assert(baumgarte_partial_dq.cols() == dimv_);
  assert(baumgarte_partial_dv.cols() == dimv_);
  assert(baumgarte_partial_da.cols() == dimv_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].computeBaumgarteDerivatives(
          model_, data_, 
          baumgarte_partial_dq.block(3*num_active_contacts, 0, 3, dimv_),
          baumgarte_partial_dv.block(3*num_active_contacts, 0, 3, dimv_),
          baumgarte_partial_da.block(3*num_active_contacts, 0, 3, dimv_));
      ++num_active_contacts;
    }
  }
}


void Robot::computeBaumgarteDerivatives(
    const double coeff, Eigen::Ref<Eigen::MatrixXd> baumgarte_partial_dq, 
    Eigen::Ref<Eigen::MatrixXd> baumgarte_partial_dv, 
    Eigen::Ref<Eigen::MatrixXd> baumgarte_partial_da) {
  assert(baumgarte_partial_dq.cols() == dimv_);
  assert(baumgarte_partial_dv.cols() == dimv_);
  assert(baumgarte_partial_da.cols() == dimv_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].computeBaumgarteDerivatives(
          model_, data_, coeff,
          baumgarte_partial_dq.block(3*num_active_contacts, 0, 3, dimv_),
          baumgarte_partial_dv.block(3*num_active_contacts, 0, 3, dimv_),
          baumgarte_partial_da.block(3*num_active_contacts, 0, 3, dimv_));
      ++num_active_contacts;
    }
  }
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


void Robot::setContactForces(const Eigen::Ref<const Eigen::VectorXd>& f) {
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].computeJointForceFromContactForce(
          f.segment<3>(3*num_active_contacts), fjoint_);
      ++num_active_contacts;
    }
    else {
      point_contacts_[i].computeJointForceFromContactForce(
          Eigen::Vector3d::Zero(), fjoint_);
    }
  }
}


void Robot::RNEA(const Eigen::Ref<const Eigen::VectorXd>& q, 
                 const Eigen::Ref<const Eigen::VectorXd>& v, 
                 const Eigen::Ref<const Eigen::VectorXd>& a, 
                 Eigen::Ref<Eigen::VectorXd> tau) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(tau.size() == dimv_);
  tau = pinocchio::rnea(model_, data_, q, v, a);
  if (point_contacts_.empty()) {
    tau = pinocchio::rnea(model_, data_, q, v, a);
  }
  else {
    tau = pinocchio::rnea(model_, data_, q, v, a, fjoint_);
  }
}


void Robot::RNEADerivatives(const Eigen::Ref<const Eigen::VectorXd>& q, 
                            const Eigen::Ref<const Eigen::VectorXd>& v, 
                            const Eigen::Ref<const Eigen::VectorXd>& a,
                            Eigen::Ref<Eigen::MatrixXd> dRNEA_partial_dq, 
                            Eigen::Ref<Eigen::MatrixXd> dRNEA_partial_dv, 
                            Eigen::Ref<Eigen::MatrixXd> dRNEA_partial_da) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(dRNEA_partial_dq.cols() == dimv_);
  assert(dRNEA_partial_dq.rows() == dimv_);
  assert(dRNEA_partial_dv.cols() == dimv_);
  assert(dRNEA_partial_dv.rows() == dimv_);
  assert(dRNEA_partial_da.cols() == dimv_);
  assert(dRNEA_partial_da.rows() == dimv_);
  if (point_contacts_.empty()) {
    pinocchio::computeRNEADerivatives(model_, data_, q, v, a, dRNEA_partial_dq, 
                                      dRNEA_partial_dv, dRNEA_partial_da);
  }
  else {
    pinocchio::computeRNEADerivatives(model_, data_, q, v, a, fjoint_, 
                                      dRNEA_partial_dq, dRNEA_partial_dv, 
                                      dRNEA_partial_da);
  }
  dRNEA_partial_da.triangularView<Eigen::StrictlyLower>() 
      = dRNEA_partial_da.transpose().triangularView<Eigen::StrictlyLower>();
}


void Robot::dRNEAPartialdFext(Eigen::Ref<Eigen::MatrixXd> dRNEA_partial_dfext) {
  assert(dRNEA_partial_dfext.rows() == dimv_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].getContactJacobian(
          model_, data_,  -1, 
          dRNEA_partial_dfext.block(0, 3*num_active_contacts, dimv_, 3), true);
      ++num_active_contacts;
    }
  }
}


void Robot::stateEquation(const Eigen::Ref<const Eigen::VectorXd>& q, 
                          const Eigen::Ref<const Eigen::VectorXd>& v, 
                          const Eigen::Ref<const Eigen::VectorXd>& tau, 
                          Eigen::Ref<Eigen::VectorXd> dq,
                          Eigen::Ref<Eigen::VectorXd> dv) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(tau.size() == dimv_);
  assert(dq.size() == dimv_);
  assert(dv.size() == dimv_);
  dq = v;
  if (point_contacts_.empty()) {
    dv = pinocchio::aba(model_, data_, q, v, tau);
  }
  else {
    dv = pinocchio::aba(model_, data_, q, v, tau, fjoint_);
  }
}


void Robot::generateFeasibleConfiguration(Eigen::Ref<Eigen::VectorXd> q) const {
  assert(q.size() == dimq_);
  Eigen::VectorXd q_min = model_.lowerPositionLimit;
  Eigen::VectorXd q_max = model_.upperPositionLimit;
  if (floating_base_.has_floating_base()) {
    q_min.head(6).array() = -1;
    q_max.head(6).array() = 1;
  }
  q = pinocchio::randomConfiguration(model_, q_min, q_max);
}


void Robot::normalizeConfiguration(Eigen::Ref<Eigen::VectorXd> q) const {
  assert(q.size() == dimq_);
  if (floating_base_.has_floating_base()) {
    if (q.segment<4>(3).squaredNorm() 
          <= std::numeric_limits<double>::epsilon()) {
      q.coeffRef(6) = 1;
    }
    pinocchio::normalize(model_, q);
  }
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