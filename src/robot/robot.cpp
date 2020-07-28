#include "robot/robot.hpp"

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
    is_each_contact_active_(),
    joint_damping_coeff_(),
    is_effective_joint_damping_(false) {
  pinocchio::urdf::buildModel(urdf_file_name, model_);
  data_ = pinocchio::Data(model_);
  floating_base_ = FloatingBase(model_);
  fjoint_ = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  dimq_ = model_.nq;
  dimv_ = model_.nv;
}


Robot::Robot(const std::string& urdf_file_name, 
             const std::vector<int>& contact_frames, 
             const double baumgarte_weight_on_position, 
             const double baumgarte_weight_on_velocity) 
  : model_(),
    data_(model_),
    urdf_file_name_(urdf_file_name),
    point_contacts_(),
    floating_base_(),
    fjoint_(),
    dimq_(0),
    dimv_(0),
    dimf_(0),
    is_each_contact_active_(),
    joint_damping_coeff_(),
    is_effective_joint_damping_(false) {
  assert(baumgarte_weight_on_position >= 0);
  assert(baumgarte_weight_on_velocity >= 0);
  pinocchio::urdf::buildModel(urdf_file_name, model_);
  data_ = pinocchio::Data(model_);
  for (const auto& frame : contact_frames) {
    point_contacts_.push_back(PointContact(model_, frame, 
                                           baumgarte_weight_on_position, 
                                           baumgarte_weight_on_velocity));
    is_each_contact_active_.push_back(false);
  }
  floating_base_ = FloatingBase(model_);
  fjoint_ = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  dimq_ = model_.nq;
  dimv_ = model_.nv;
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
    is_each_contact_active_(),
    joint_damping_coeff_(),
    is_effective_joint_damping_(false) {
}


Robot::~Robot() {
}


void Robot::buildRobotModelFromXML(const std::string& xml) {
  pinocchio::urdf::buildModelFromXML(xml, model_);
  data_ = pinocchio::Data(model_);
  floating_base_ = FloatingBase(model_);
  fjoint_ = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  dimq_ = model_.nq;
  dimv_ = model_.nv;
}


void Robot::buildRobotModelFromXML(const std::string& xml,
                                   const std::vector<int>& contact_frames, 
                                   const double baumgarte_weight_on_position, 
                                   const double baumgarte_weight_on_velocity) {
  pinocchio::urdf::buildModelFromXML(xml, model_);
  data_ = pinocchio::Data(model_);
  point_contacts_.clear();
  is_each_contact_active_.clear();
  for (const auto& frame : contact_frames) {
    point_contacts_.push_back(PointContact(model_, frame, 
                                           baumgarte_weight_on_position, 
                                           baumgarte_weight_on_velocity));
    is_each_contact_active_.push_back(false);
  }
  floating_base_ = FloatingBase(model_);
  fjoint_ = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  dimq_ = model_.nq;
  dimv_ = model_.nv;
}


void Robot::integrateConfiguration(const Eigen::VectorXd& v, 
                                   const double integration_length, 
                                   Eigen::VectorXd& q) const {
  assert(v.size() == dimv_);
  assert(integration_length >= 0);
  assert(q.size() == dimq_);
  q = pinocchio::integrate(model_, q, integration_length*v);
}


void Robot::differenceConfiguration(const Eigen::VectorXd& q_plus, 
                                    const Eigen::VectorXd& q_minus,
                                    Eigen::VectorXd& difference) const {
  assert(q_plus.size() == dimq_);
  assert(q_minus.size() == dimq_);
  assert(difference.size() == dimv_);
  pinocchio::difference(model_, q_minus, q_plus, difference);
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
  pinocchio::dIntegrate(model_, q, integration_length*v, dIntegrate_dq, 
                        pinocchio::ARG0);
  pinocchio::dIntegrate(model_, q, integration_length*v, dIntegrate_dv, 
                        pinocchio::ARG1);
}


void Robot::updateKinematics(const Eigen::VectorXd& q, const Eigen::VectorXd& v,
                             const Eigen::VectorXd& a) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  pinocchio::forwardKinematics(model_, data_, q, v, a);
  pinocchio::updateFramePlacements(model_, data_);
  pinocchio::computeForwardKinematicsDerivatives(model_, data_, q, v, a);
  for (int i=0; i<point_contacts_.size(); ++i) {
    point_contacts_[i].resetContactPointByCurrentKinematics(data_);
  }
}


void Robot::computeBaumgarteResidual(
    Eigen::VectorXd& baumgarte_residual) const {
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].computeBaumgarteResidual(model_, data_, 
                                                  3*num_active_contacts, 
                                                  baumgarte_residual);
      ++num_active_contacts;
    }
  }
  assert(baumgarte_residual.size() >= 3*num_active_contacts);
}


void Robot::computeBaumgarteDerivatives(Eigen::MatrixXd& baumgarte_partial_dq, 
                                        Eigen::MatrixXd& baumgarte_partial_dv,
                                        Eigen::MatrixXd& baumgarte_partial_da) {
  assert(baumgarte_partial_dq.cols() == dimv_);
  assert(baumgarte_partial_dv.cols() == dimv_);
  assert(baumgarte_partial_da.cols() == dimv_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].computeBaumgarteDerivatives(model_, data_, 
                                                     3*num_active_contacts, 
                                                     baumgarte_partial_dq, 
                                                     baumgarte_partial_dv, 
                                                     baumgarte_partial_da);
      ++num_active_contacts;
    }
  }
  assert(baumgarte_partial_dq.rows() >= 3*num_active_contacts);
  assert(baumgarte_partial_dv.rows() >= 3*num_active_contacts);
  assert(baumgarte_partial_da.rows() >= 3*num_active_contacts);
}


void Robot::setActiveContacts(const std::vector<bool>& is_each_contact_active) {
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
  dimf_ = 3*num_active_contacts;
}


void Robot::setContactForces(const Eigen::VectorXd& fext) {
  assert(fext.size() <= 3*point_contacts_.size());
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
  if (is_effective_joint_damping_) {
    tau.array() += joint_damping_coeff_.array() * v.array();
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
  if (is_effective_joint_damping_) {
    for (int i=0; i<dimv_; ++i) {
      dRNEA_partial_dv.coeffRef(i, i) += joint_damping_coeff_.coeff(i);
    }
  }
}


void Robot::dRNEAPartialdFext(Eigen::MatrixXd& dRNEA_partial_dfext) {
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].getContactJacobian(model_, data_,  
                                            3*num_active_contacts, 
                                            dRNEA_partial_dfext);
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
  if (is_effective_joint_damping_) {
    u.array() -= joint_damping_coeff_.array() * v.array();
  }
  dv = pinocchio::aba(model_, data_, q, v, u);
}


void Robot::setPassiveTorques(Eigen::VectorXd& torques) const {
  floating_base_.setPassiveTorques(torques);
}


void Robot::passiveConstraintViolation(const Eigen::VectorXd& torques, 
                                       Eigen::VectorXd& violation) const {
  floating_base_.computePassiveConstraintViolation(torques, violation);
}


bool Robot::hasFloatingBase() {
  return floating_base_.has_floating_base();
}


Eigen::VectorXd Robot::jointEffortLimit() const {
  return model_.effortLimit;
}


Eigen::VectorXd Robot::jointVelocityLimit() const {
  return model_.velocityLimit;
}


Eigen::VectorXd Robot::lowerJointPositionLimit() const {
  return model_.lowerPositionLimit;
}


Eigen::VectorXd Robot::upperJointPositionLimit() const {
  return model_.upperPositionLimit;
}


void Robot::setJointDamping(const Eigen::VectorXd& joint_damping_coeff) {
  assert(joint_damping_coeff.size() == dimv_);
  assert(joint_damping_coeff.minCoeff() >= 0);
  joint_damping_coeff_ = joint_damping_coeff;
  is_effective_joint_damping_ = true;
}


int Robot::dimq() const {
  return dimq_;
}


int Robot::dimv() const {
  return dimv_;
}


int Robot::dim_passive() const {
  return floating_base_.dim_passive();
}


int Robot::dimf() const {
  return dimf_;
}


int Robot::max_dimf() const {
  return 3*point_contacts_.size();
}


int Robot::max_point_contacts() const {
  return point_contacts_.size();
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