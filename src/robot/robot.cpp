#include "robot/robot.hpp"

#include <algorithm>


namespace invdynocp {

Robot::Robot(const std::string& urdf_file_name, 
             const unsigned int max_point_contacts)
  : model_(),
    data_(model_),
    urdf_file_name_(urdf_file_name),
    point_contacts_(),
    passive_joints_(),
    fjoint_(),
    dimq_(0),
    dimv_(0),
    max_point_contacts_(max_point_contacts) {
  // Build Pinocchio model from URDF.
  pinocchio::urdf::buildModel(urdf_file_name, model_);
  data_ = pinocchio::Data(model_);
  passive_joints_ = PassiveJoints(model_);
  fjoint_ = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  dimq_ = model_.nq;
  dimv_ = model_.nv;
}


Robot::~Robot() {
}


Robot::Robot(const Robot& other) 
  : model_(), 
    data_(model_),
    urdf_file_name_(other.urdf_file_name()),
    point_contacts_(),
    passive_joints_(),
    fjoint_(),
    dimq_(0),
    dimv_(0),
    max_point_contacts_(other.max_point_contacts()) {
  // Build Pinocchio model from URDF.
  pinocchio::urdf::buildModel(other.urdf_file_name(), model_);
  data_ = pinocchio::Data(model_);
  passive_joints_ = PassiveJoints(model_);
  fjoint_ = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  dimq_ = model_.nq;
  dimv_ = model_.nv;
}


Robot::Robot& operator=(const Robot& other) {
  // Build Pinocchio model from URDF.
  pinocchio::urdf::buildModel(other.urdf_file_name(), model_);
  data_ = pinocchio::Data(model_);
  urdf_file_name_ = other.urdf_file_name();
  passive_joints_ = PassiveJoints(model_);
  fjoint_ = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  dimq_ = model_.nq;
  dimv_ = model_.nv;
  max_point_contacts_ = other.max_point_contacts();
}


Robot::Robot(const Robot&& other) noexcept
  : model_(), 
    data_(model_),
    urdf_file_name_(other.urdf_file_name()),
    point_contacts_(),
    passive_joints_(),
    fjoint_(),
    dimq_(0),
    dimv_(0),
    max_point_contacts_(other.max_point_contacts()) {
  // Build Pinocchio model from URDF.
  pinocchio::urdf::buildModel(other.urdf_file_name(), model_);
  data_ = pinocchio::Data(model_);
  passive_joints_ = PassiveJoints(model_);
  fjoint_ = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  dimq_ = model_.nq;
  dimv_ = model_.nv;
}


Robot::Robot& operator=(const Robot&& other) noexcept {
  // Build Pinocchio model from URDF.
  pinocchio::urdf::buildModel(other.urdf_file_name(), model_);
  data_ = pinocchio::Data(model_);
  urdf_file_name_ = other.urdf_file_name();
  passive_joints_ = PassiveJoints(model_);
  fjoint_ = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  dimq_ = model_.nq;
  dimv_ = model_.nv;
  max_point_contacts_ = other.max_point_contacts();
}


void Robot::integrateConfiguration(const Eigen::VectorXd& v, 
                                   const double integration_length,
                                   Eigen::VectorXd& q) {
  q = pinocchio::integrate(model_, q, integration_length*v);
}


void Robot::integrateConfiguration(const Eigen::VectorXd& q, 
                                   const Eigen::VectorXd& v, 
                                   const double integration_length, 
                                   Eigen::VectorXd& q_plus) {
  pinocchio::integrate(model_, q, integration_length*v, q_plus);
}


void Robot::differenceConfigurations(const Eigen::VectorXd& q_plus, 
                                     const Eigen::VectorXd& q_minus, 
                                     Eigen::VectorXd& difference) {
  difference = pinocchio::difference(model_, q_minus, q_plus);
}


void Robot::updateKinematics(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                             const Eigen::VectorXd& a) {
  pinocchio::forwardKinematics(model_, data_, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model_, data_, q, v, a);
  pinocchio::updateFramePlacements(model_, data_);
  for (int i=0; i<contacts_.size(); ++i) {
    contacts_[i].resetContactPointByCurrentKinematics(data_);
  }
}


void Robot::computeBaumgarteResidual(Eigen::VectorXd& baumgarte_residual) {
  for (int i=0; i<contacts_.size(); ++i) {
    contacts_[i].computeBaumgarteResidual(model_, data_, 3*i, 
                                          baumgarte_residual);
  }
}


void Robot::computeBaumgarteDerivatives(Eigen::MatrixXd& baumgarte_partial_dq, 
                                        Eigen::MatrixXd& baumgarte_partial_dv, 
                                        Eigen::MatrixXd& baumgarte_partial_da) {
  for (int i=0; i<contacts_.size(); ++i) {
    contacts_[i].computeBaumgarteDerivatives(model_, data_, 0, 3*i, 
                                             dBaumgarte_dq_, dBaumgarte_dv_, 
                                             dBaumgarte_da_);
  }
}


void Robot::setFext(const Eigen::VectorXd& fext) {
  for (int i=0; i<contacts_.size(); ++i) {
    contacts_[i].computeJointForceFromContactForce(fext.segment<3>(3*i), 
                                                   fjoint_);
  }
}


void Robot::RNEA(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                 const Eigen::VectorXd& a, Eigen::VectorXd& tau) {
  tau = pinocchio::rnea(model_, data_, q, v, a);
}


void Robot::RNEA(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                 const Eigen::VectorXd& a, const Eigen::VectorXd& fext, 
                 Eigen::VectorXd& tau) {
  tau = pinocchio::rnea(model_, data_, q, v, a, fjoint_);
}


void Robot::RNEADerivatives(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                            const Eigen::VectorXd& a, 
                            Eigen::MatrixXd& dRNEA_partial_dq, 
                            Eigen::MatrixXd& dRNEA_partial_dv, 
                            Eigen::MatrixXd& dRNEA_partial_da) {
  pinocchio::computeRNEADerivatives(model_, data_, q, v, a, dRNEA_partial_dq,
                                    dRNEA_partial_dv, dRNEA_partial_da);
  data_.M.triangularView<Eigen::StrictlyLower>() 
      = data_.M.transpose().triangularView<Eigen::StrictlyLower>();
}


void Robot::RNEADerivatives(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                            const Eigen::VectorXd& a, 
                            const Eigen::VectorXd& fext, 
                            Eigen::MatrixXd& dRNEA_partial_dq, 
                            Eigen::MatrixXd& dRNEA_partial_dv, 
                            Eigen::MatrixXd& dRNEA_partial_da_and_fext) {
  pinocchio::computeRNEADerivatives(model_, data_, q, v, a, dRNEA_partial_dq,
                                    dRNEA_partial_dv, 
                                    dRNEA_partial_da_and_fext.topRows(dimv_));
  data_.M.triangularView<Eigen::StrictlyLower>() 
      = data_.M.transpose().triangularView<Eigen::StrictlyLower>();
  for (int i=0; i<contacts_.size(); ++i) {
    contacts_[i].computeContactJacobian(model_, data_, 0, dimv_+3*i, 
                                        dRNEA_partial_da_and_fext);
  }
}


void Robot::addPointContact(const unsigned int contact_frame_id, 
                            const double baumgarte_alpha, 
                            const double baumgarte_beta) {
  auto find
      = std::find_if(contacts_.begin(), contacts_.end(), 
                     [contact_frame_id](PointContact& contact)
                     { return contact.contact_frame_id()==contact_frame_id; });
  if (find == contacts_.end()) {
    contacts_.push_back(PointContact(model_, contact_frame_id, baumgarte_alpha, 
                                     baumgarte_beta));
  }
}


void Robot::removePointContact(const unsigned int contact_frame_id) {
  auto find
      = std::find_if(contacts_.begin(), contacts_.end(), 
                     [contact_frame_id](PointContact& contact)
                     { return contact.contact_frame_id()==contact_frame_id; });
  if (find != contacts_.end()) {
    std::swap<PointContact>(*find, contacts_.back());
    contacts_.pop_back();
  }
}


void Robot::setPassiveTorques(Eigen::VectorXd& tau) const {
  passive_joints_.setPassiveTorques(tau);
}


void Robot::passiveConstraintViolation(const Eigen::VectorXd& tau, 
                                       Eigen::VectorXd& violation) const {
  passive_joints_.computePassiveConstraintViolation(tau, residual);
}


void Robot::passiveConstraintsDerivative(Eigen::MatrixXd& derivative) const {
  passive_joints_.computePassiveConstraintDerivative(derivative);
}


unsigned int Robot::dimq() const {
  return dimq_;
}


unsigned int Robot::dimv() const {
  return dimv_;
}


unsigned int Robot::dimf() const {
  return 3*contacts_.size();
}


unsigned int Robot::dimfmax() const {
  return 3*max_point_contacts_;
}


unsigned int Robot::dim_passive() const {
  return passive_joints_.dim_passive();
}


unsigned int Robot::max_point_contacts() const {
  return max_point_contacts_;
}

} // namespace invdynocp