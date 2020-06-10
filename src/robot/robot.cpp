#include "robot/robot.hpp"

#include <assert.h>


namespace idocp {

Robot::Robot(const std::string& urdf_file_name)
  : model_(),
    data_(model_),
    urdf_file_name_(urdf_file_name),
    point_contacts_(),
    passive_joints_(),
    fjoint_(),
    dimq_(0),
    dimv_(0) {
  // Build Pinocchio model from URDF.
  pinocchio::urdf::buildModel(urdf_file_name, model_);
  data_ = pinocchio::Data(model_);
  passive_joints_ = PassiveJoints(model_);
  fjoint_ = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  dimq_ = model_.nq;
  dimv_ = model_.nv;
}


Robot::Robot(const std::string& urdf_file_name, 
             const std::vector<unsigned int>& contact_frames, 
             const double baumgarte_weight_on_position, 
             const double baumgarte_weight_on_velocity) 
  : model_(),
    data_(model_),
    urdf_file_name_(urdf_file_name),
    point_contacts_(),
    passive_joints_(),
    fjoint_(),
    dimq_(0),
    dimv_(0) {
  assert(baumgarte_weight_on_position >= 0);
  assert(baumgarte_weight_on_velocity >= 0);
  // Build Pinocchio model from URDF.
  pinocchio::urdf::buildModel(urdf_file_name, model_);
  data_ = pinocchio::Data(model_);
  for (const auto& frame : contact_frames) {
    point_contacts_.push_back(PointContact(model_, frame, 
                                           baumgarte_weight_on_position, 
                                           baumgarte_weight_on_velocity));
  }
  passive_joints_ = PassiveJoints(model_);
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
    passive_joints_(),
    fjoint_(),
    dimq_(0),
    dimv_(0) {
}


Robot::~Robot() {
}


void Robot::integrateConfiguration(const Eigen::VectorXd& v, 
                                   const double integration_length, 
                                   Eigen::VectorXd& q) const {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  q =  pinocchio::integrate(model_, q, integration_length*v);
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
  pinocchio::dIntegrate(model_, q, v, dIntegrate_dq, pinocchio::ARG0);
  pinocchio::dIntegrate(model_, q, v, dIntegrate_dv, pinocchio::ARG1);
}


void Robot::updateKinematics(const Eigen::VectorXd& q, const Eigen::VectorXd& v,
                             const Eigen::VectorXd& a) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  pinocchio::forwardKinematics(model_, data_, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model_, data_, q, v, a);
  pinocchio::updateFramePlacements(model_, data_);
  for (int i=0; i<point_contacts_.size(); ++i) {
    point_contacts_[i].resetContactPointByCurrentKinematics(data_);
  }
}


void Robot::computeBaumgarteResidual(
    Eigen::VectorXd& baumgarte_residual) const {
  unsigned int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].computeBaumgarteResidual(model_, data_, 
                                                  3*num_active_contacts, 
                                                  baumgarte_residual);
      ++num_active_contacts;
    }
  }
  assert(3*num_active_contacts == baumgarte_residual.size());
}


void Robot::computeBaumgarteDerivatives(Eigen::MatrixXd& baumgarte_partial_dq, 
                                        Eigen::MatrixXd& baumgarte_partial_dv,
                                        Eigen::MatrixXd& baumgarte_partial_da) {
  unsigned int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].computeBaumgarteDerivatives(model_, data_, 0, 
                                                     3*num_active_contacts, 
                                                     baumgarte_partial_dq, 
                                                     baumgarte_partial_dv, 
                                                     baumgarte_partial_da);
      ++num_active_contacts;
    }
  }
}


void Robot::setActiveContacts(const std::vector<bool>& is_each_contact_active, 
                              const Eigen::VectorXd& fext) {
  assert(is_each_contact_active.size() == point_contacts_.size());
  unsigned int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if(is_each_contact_active[i]) {
      point_contacts_[i].activate();
      point_contacts_[i].computeJointForceFromContactForce(
          fext.segment<3>(3*num_active_contacts), fjoint_);
      ++num_active_contacts;
    }
    else {
      point_contacts_[i].deactivate();
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
  unsigned int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].getContactJacobian(model_, data_,  
                                            3*num_active_contacts, 0, 
                                            dRNEA_partial_dfext);
      ++num_active_contacts;
    }
  }
}


void Robot::setPassiveTorques(Eigen::VectorXd& tau) const {
  passive_joints_.setPassiveTorques(tau);
}


void Robot::passiveConstraintViolation(const Eigen::VectorXd& tau, 
                                       Eigen::VectorXd& violation) const {
  passive_joints_.computePassiveConstraintViolation(tau, violation);
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


unsigned int Robot::dimq() const {
  return dimq_;
}


unsigned int Robot::dimv() const {
  return dimv_;
}


unsigned int Robot::dim_passive() const {
  return passive_joints_.dim_passive();
}


unsigned int Robot::max_point_contacts() const {
  return point_contacts_.size();
}

} // namespace idocp 