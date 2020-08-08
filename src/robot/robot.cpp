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
    dimf_(0),
    max_dimf_(0),
    is_each_contact_active_(),
    joint_effort_limit_(),
    joint_velocity_limit_(),
    lower_joint_position_limit_(),
    upper_joint_position_limit_(),
    configuration_jacobian_() {
  pinocchio::urdf::buildModel(urdf_file_name, model_);
  data_ = pinocchio::Data(model_);
  fjoint_ = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  floating_base_ = FloatingBase(model_);
  dimq_ = model_.nq;
  dimv_ = model_.nv;
  const int dim_joint = model_.nv - floating_base_.dim_passive();
  joint_effort_limit_.resize(dim_joint);
  joint_velocity_limit_.resize(dim_joint);
  lower_joint_position_limit_.resize(dim_joint);
  upper_joint_position_limit_.resize(dim_joint);
  joint_effort_limit_ = model_.effortLimit.tail(dim_joint);
  joint_velocity_limit_ = model_.velocityLimit.tail(dim_joint);
  lower_joint_position_limit_ = model_.lowerPositionLimit.tail(dim_joint);
  upper_joint_position_limit_ = model_.upperPositionLimit.tail(dim_joint);
  configuration_jacobian_.resize(model_.nq, model_.nv);
  configuration_jacobian_.fill(0);
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
    dimf_(0),
    max_dimf_(0),
    is_each_contact_active_(),
    joint_effort_limit_(),
    joint_velocity_limit_(),
    lower_joint_position_limit_(),
    upper_joint_position_limit_(),
    configuration_jacobian_() {
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
  const int dim_joint = model_.nv - floating_base_.dim_passive();
  joint_effort_limit_.resize(dim_joint);
  joint_velocity_limit_.resize(dim_joint);
  lower_joint_position_limit_.resize(dim_joint);
  upper_joint_position_limit_.resize(dim_joint);
  joint_effort_limit_ = model_.effortLimit.tail(dim_joint);
  joint_velocity_limit_ = model_.velocityLimit.tail(dim_joint);
  lower_joint_position_limit_ = model_.lowerPositionLimit.tail(dim_joint);
  upper_joint_position_limit_ = model_.upperPositionLimit.tail(dim_joint);
  configuration_jacobian_.resize(model_.nq, model_.nv);
  configuration_jacobian_.fill(0);
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
    dimf_(0),
    max_dimf_(0),
    is_each_contact_active_(),
    joint_effort_limit_(),
    joint_velocity_limit_(),
    lower_joint_position_limit_(),
    upper_joint_position_limit_(),
    configuration_jacobian_() {
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
  const int dim_joint = model_.nv - floating_base_.dim_passive();
  joint_effort_limit_.resize(dim_joint);
  joint_velocity_limit_.resize(dim_joint);
  lower_joint_position_limit_.resize(dim_joint);
  upper_joint_position_limit_.resize(dim_joint);
  joint_effort_limit_ = model_.effortLimit.tail(dim_joint);
  joint_velocity_limit_ = model_.velocityLimit.tail(dim_joint);
  lower_joint_position_limit_ = model_.lowerPositionLimit.tail(dim_joint);
  upper_joint_position_limit_ = model_.upperPositionLimit.tail(dim_joint);
  configuration_jacobian_.resize(model_.nq, model_.nv);
  configuration_jacobian_.fill(0);
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
  const int dim_joint = model_.nv - floating_base_.dim_passive();
  joint_effort_limit_.resize(dim_joint);
  joint_velocity_limit_.resize(dim_joint);
  lower_joint_position_limit_.resize(dim_joint);
  upper_joint_position_limit_.resize(dim_joint);
  joint_effort_limit_ = model_.effortLimit.tail(dim_joint);
  joint_velocity_limit_ = model_.velocityLimit.tail(dim_joint);
  lower_joint_position_limit_ = model_.lowerPositionLimit.tail(dim_joint);
  upper_joint_position_limit_ = model_.upperPositionLimit.tail(dim_joint);
  configuration_jacobian_.resize(model_.nq, model_.nv);
  configuration_jacobian_.fill(0);
}


void Robot::integrateConfiguration(const Eigen::VectorXd& v, 
                                   const double integration_length, 
                                   Eigen::VectorXd& q) const {
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


void Robot::subtractConfiguration(const Eigen::VectorXd& q_plus, 
                                  const Eigen::VectorXd& q_minus,
                                  Eigen::VectorXd& difference) const {
  assert(q_plus.size() == dimq_);
  assert(q_minus.size() == dimq_);
  assert(difference.size() == dimv_);
  if (floating_base_.has_floating_base()) {
    pinocchio::difference(model_, q_minus, q_plus, difference);
  }
  else {
    difference = q_plus - q_minus;
  }
}


void Robot::dIntegrateConfiguration(const Eigen::VectorXd& q, 
                                    const Eigen::VectorXd& v,
                                    const double integration_length,
                                    Eigen::MatrixXd& dIntegrate_dq,
                                    Eigen::MatrixXd& dIntegrate_dv) const {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(integration_length >= 0);
  assert(dIntegrate_dq.rows() == dimv_);
  assert(dIntegrate_dq.cols() == dimv_);
  assert(dIntegrate_dv.rows() == dimv_);
  assert(dIntegrate_dv.cols() == dimv_);
  pinocchio::dIntegrate(model_, q, v, dIntegrate_dq, 
                        pinocchio::ARG0);
  pinocchio::dIntegrate(model_, q, v, dIntegrate_dv, 
                        pinocchio::ARG1);
}


void Robot::computeConfigurationJacobian(const Eigen::VectorXd& q) {
  assert(q.size() == dimq_);
  pinocchio::integrateCoeffWiseJacobian(model_, q, configuration_jacobian_);
}


void Robot::computeTangentGradient(
    const Eigen::VectorXd& gradient_at_configuration, 
    Eigen::VectorXd& gradient_at_tangent) const {
  assert(gradient_at_configuration.size() == dimq_);
  assert(gradient_at_tangent.size() == dimv_);
  gradient_at_tangent 
      = configuration_jacobian_.transpose() * gradient_at_configuration;
}


void Robot::computeTangentHessian(
    const Eigen::MatrixXd& hessian_at_configuration, 
    Eigen::MatrixXd& hessian_at_tangent) const {
  assert(hessian_at_configuration.rows() == dimq_);
  assert(hessian_at_configuration.cols() == dimq_);
  assert(hessian_at_tangent.rows() == dimv_);
  assert(hessian_at_tangent.cols() == dimv_);
  hessian_at_tangent 
      = configuration_jacobian_.transpose() * hessian_at_configuration
                                            * configuration_jacobian_;
}


void Robot::augmentTangentHessian(
    const Eigen::MatrixXd& hessian_at_configuration, const double coeff,
    Eigen::MatrixXd& augmented_hessian_at_tangent) const {
  assert(hessian_at_configuration.rows() == dimq_);
  assert(hessian_at_configuration.cols() == dimq_);
  assert(augmented_hessian_at_tangent.rows() == dimv_);
  assert(augmented_hessian_at_tangent.cols() == dimv_);
  augmented_hessian_at_tangent.noalias()
      += coeff * configuration_jacobian_.transpose() * hessian_at_configuration
                                                     * configuration_jacobian_;
}


void Robot::updateKinematics(const Eigen::VectorXd& q, const Eigen::VectorXd& v,
                             const Eigen::VectorXd& a) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  pinocchio::forwardKinematics(model_, data_, q, v, a);
  pinocchio::updateFramePlacements(model_, data_);
  pinocchio::computeForwardKinematicsDerivatives(model_, data_, q, v, a);
}


void Robot::computeBaumgarteResidual(const int block_begin, 
                                     Eigen::VectorXd& baumgarte_residual) const {
  assert(baumgarte_residual.size() >= max_dimf_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].computeBaumgarteResidual(
          model_, data_, block_begin+3*num_active_contacts, baumgarte_residual);
      ++num_active_contacts;
    }
  }
}


void Robot::computeBaumgarteResidual(const int block_begin, const double coeff, 
                                     Eigen::VectorXd& baumgarte_residual) const {
  assert(baumgarte_residual.size() >= max_dimf_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].computeBaumgarteResidual(
          model_, data_, block_begin+3*num_active_contacts, coeff, 
          baumgarte_residual);
      ++num_active_contacts;
    }
  }
}


void Robot::computeBaumgarteDerivatives(const int block_rows_begin,
                                        Eigen::MatrixXd& baumgarte_partial_dq, 
                                        Eigen::MatrixXd& baumgarte_partial_dv,
                                        Eigen::MatrixXd& baumgarte_partial_da) {
  assert(baumgarte_partial_dq.cols() == dimv_);
  assert(baumgarte_partial_dq.rows() == block_rows_begin+max_dimf_);
  assert(baumgarte_partial_dv.cols() == dimv_);
  assert(baumgarte_partial_dv.rows() == block_rows_begin+max_dimf_);
  assert(baumgarte_partial_da.cols() == dimv_);
  assert(baumgarte_partial_da.rows() == block_rows_begin+max_dimf_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].computeBaumgarteDerivatives(
          model_, data_, block_rows_begin+3*num_active_contacts, 
          baumgarte_partial_dq, baumgarte_partial_dv, baumgarte_partial_da);
      ++num_active_contacts;
    }
  }
}


void Robot::computeBaumgarteDerivatives(const int block_rows_begin, 
                                        const double coeff,
                                        Eigen::MatrixXd& baumgarte_partial_dq, 
                                        Eigen::MatrixXd& baumgarte_partial_dv,
                                        Eigen::MatrixXd& baumgarte_partial_da) {
  assert(baumgarte_partial_dq.cols() == dimv_);
  assert(baumgarte_partial_dq.rows() == block_rows_begin+max_dimf_);
  assert(baumgarte_partial_dv.cols() == dimv_);
  assert(baumgarte_partial_dv.rows() == block_rows_begin+max_dimf_);
  assert(baumgarte_partial_da.cols() == dimv_);
  assert(baumgarte_partial_da.rows() == block_rows_begin+max_dimf_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].computeBaumgarteDerivatives(
          model_, data_, block_rows_begin+3*num_active_contacts, coeff,
          baumgarte_partial_dq, baumgarte_partial_dv, baumgarte_partial_da);
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
  dimf_ = 3 * num_active_contacts;
}


void Robot::setContactForces(const Eigen::VectorXd& fext) {
  assert(fext.size() == max_dimf_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].computeJointForceFromContactForce(
          fext.segment<3>(3*num_active_contacts), fjoint_);
      ++num_active_contacts;
    }
    else {
      point_contacts_[i].computeJointForceFromContactForce(
          Eigen::Vector3d::Zero(), fjoint_);
    }
  }
}


void Robot::RNEA(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                 const Eigen::VectorXd& a, Eigen::VectorXd& tau) {
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


void Robot::RNEADerivatives(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                            const Eigen::VectorXd& a, 
                            Eigen::MatrixXd& dRNEA_partial_dq, 
                            Eigen::MatrixXd& dRNEA_partial_dv, 
                            Eigen::MatrixXd& dRNEA_partial_da) {
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


void Robot::dRNEAPartialdFext(Eigen::MatrixXd& dRNEA_partial_dfext) {
  assert(dRNEA_partial_dfext.cols() == max_dimf_);
  assert(dRNEA_partial_dfext.rows() == dimv_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].getContactJacobian(model_, data_,  
                                            3*num_active_contacts, -1,
                                            dRNEA_partial_dfext, true);
      ++num_active_contacts;
    }
  }
}


void Robot::stateEquation(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                          const Eigen::VectorXd& tau, Eigen::VectorXd& dq, 
                          Eigen::VectorXd& dv) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(tau.size() == dimv_);
  assert(dq.size() == dimv_);
  assert(dv.size() == dimv_);
  dq = v;
  Eigen::VectorXd u = tau;
  if (point_contacts_.empty()) {
    dv = pinocchio::aba(model_, data_, q, v, u);
  }
  else {
    dv = pinocchio::aba(model_, data_, q, v, u, fjoint_);
  }
}


void Robot::generateFeasibleConfiguration(Eigen::VectorXd& q) const {
  assert(q.size() == dimq_);
  Eigen::VectorXd q_min = model_.lowerPositionLimit;
  Eigen::VectorXd q_max = model_.upperPositionLimit;
  if (floating_base_.has_floating_base()) {
    q_min.head(6).array() = -1;
    q_max.head(6).array() = 1;
  }
  q = pinocchio::randomConfiguration(model_, q_min, q_max);
}


void Robot::normalizeConfiguration(Eigen::VectorXd& q) const {
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


int Robot::dimf() const {
  return dimf_;
}


int Robot::max_dimf() const {
  return max_dimf_;
}


bool Robot::has_floating_base() const {
  return floating_base_.has_floating_base();
}


int Robot::dim_passive() const {
  return floating_base_.dim_passive();
}


std::vector<int> Robot::passive_joint_indices() const {
  return floating_base_.passive_joint_indices();
}


int Robot::max_point_contacts() const {
  return point_contacts_.size();
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

} // namespace idocp 