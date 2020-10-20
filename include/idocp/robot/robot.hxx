#ifndef IDOCP_ROBOT_HXX_
#define IDOCP_ROBOT_HXX_

#include "idocp/robot/robot.hpp"

#include <stdexcept>
#include <assert.h>

namespace idocp {

template <typename TangentVectorType, typename ConfigVectorType>
inline void Robot::integrateConfiguration(
    const Eigen::MatrixBase<TangentVectorType>& v, 
    const double integration_length, 
    const Eigen::MatrixBase<ConfigVectorType>& q) const {
  assert(v.size() == dimv_);
  assert(q.size() == dimq_);
  if (floating_base_.has_floating_base()) {
    const Eigen::VectorXd q_tmp = q;
    pinocchio::integrate(model_, q_tmp, integration_length*v, 
                         const_cast<Eigen::MatrixBase<ConfigVectorType>&>(q));
  }
  else {
    (const_cast<Eigen::MatrixBase<ConfigVectorType>&>(q)).noalias() 
        += integration_length * v;
  }
}


template <typename ConfigVectorType1, typename TangentVectorType,  
          typename ConfigVectorType2>
inline void Robot::integrateConfiguration(
    const Eigen::MatrixBase<ConfigVectorType1>& q, 
    const Eigen::MatrixBase<TangentVectorType>& v, 
    const double integration_length, 
    const Eigen::MatrixBase<ConfigVectorType2>& q_integrated) const {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(q_integrated.size() == dimq_);
  if (floating_base_.has_floating_base()) {
    pinocchio::integrate(
        model_, q, integration_length*v, 
        const_cast<Eigen::MatrixBase<ConfigVectorType2>&>(q_integrated));
  }
  else {
    (const_cast<Eigen::MatrixBase<ConfigVectorType2>&>(q_integrated)).noalias() 
        = q + integration_length * v;
  }
}


template <typename ConfigVectorType1, typename ConfigVectorType2, 
          typename TangentVectorType>
inline void Robot::subtractConfiguration(
    const Eigen::MatrixBase<ConfigVectorType1>& q_plus, 
    const Eigen::MatrixBase<ConfigVectorType2>& q_minus,
    const Eigen::MatrixBase<TangentVectorType>& difference) const {
  assert(q_plus.size() == dimq_);
  assert(q_minus.size() == dimq_);
  assert(difference.size() == dimv_);
  if (floating_base_.has_floating_base()) {
    pinocchio::difference(
        model_, q_minus, q_plus, 
        const_cast<Eigen::MatrixBase<TangentVectorType>&>(difference));
  }
  else {
    const_cast<Eigen::MatrixBase<TangentVectorType>&>(difference) 
        = q_plus - q_minus;
  }
}


template <typename ConfigVectorType1, typename ConfigVectorType2, 
          typename MatrixType>
inline void Robot::dSubtractdConfigurationPlus(
    const Eigen::MatrixBase<ConfigVectorType1>& q_plus,
    const Eigen::MatrixBase<ConfigVectorType2>& q_minus,
    const Eigen::MatrixBase<MatrixType>& dSubtract_dqplus) const {
  assert(q_plus.size() == dimq_);
  assert(q_minus.size() == dimq_);
  assert(dSubtract_dqplus.rows() == dimv_);
  assert(dSubtract_dqplus.cols() == dimv_);
  pinocchio::dDifference(
      model_, q_minus, q_plus, 
      const_cast<Eigen::MatrixBase<MatrixType>&>(dSubtract_dqplus),
      pinocchio::ARG1);
}


template <typename ConfigVectorType1, typename ConfigVectorType2, 
          typename MatrixType>
inline void Robot::dSubtractdConfigurationMinus(
    const Eigen::MatrixBase<ConfigVectorType1>& q_plus,
    const Eigen::MatrixBase<ConfigVectorType2>& q_minus,
    const Eigen::MatrixBase<MatrixType>& dSubtract_dqminus) const {
  assert(q_plus.size() == dimq_);
  assert(q_minus.size() == dimq_);
  assert(dSubtract_dqminus.rows() == dimv_);
  assert(dSubtract_dqminus.cols() == dimv_);
  pinocchio::dDifference(
      model_, q_minus, q_plus, 
      const_cast<Eigen::MatrixBase<MatrixType>&>(dSubtract_dqminus),
      pinocchio::ARG0);
}


inline const Eigen::Vector3d& Robot::framePosition(const int frame_id) const {
  return data_.oMf[frame_id].translation();
}


inline const Eigen::Matrix3d& Robot::frameRotation(const int frame_id) const {
  return data_.oMf[frame_id].rotation();
}


inline const pinocchio::SE3& Robot::framePlacement(const int frame_id) const {
  return data_.oMf[frame_id];
}


template <typename MatrixType>
inline void Robot::getFrameJacobian(const int frame_id, 
                                    const Eigen::MatrixBase<MatrixType>& J) {
  assert(J.rows() == 6);
  assert(J.cols() == dimv_);
  pinocchio::getFrameJacobian(model_, data_, frame_id, pinocchio::LOCAL, 
                              const_cast<Eigen::MatrixBase<MatrixType>&>(J));
}


template <typename ConfigVectorType, typename TangentVectorType1, 
          typename TangentVectorType2>
inline void Robot::updateKinematics(
    const Eigen::MatrixBase<ConfigVectorType>& q, 
    const Eigen::MatrixBase<TangentVectorType1>& v, 
    const Eigen::MatrixBase<TangentVectorType2>& a) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  pinocchio::forwardKinematics(model_, data_, q, v, a);
  pinocchio::updateFramePlacements(model_, data_);
  pinocchio::computeForwardKinematicsDerivatives(model_, data_, q, v, a);
}


template <typename ConfigVectorType, typename TangentVectorType>
inline void Robot::updateKinematics(
    const Eigen::MatrixBase<ConfigVectorType>& q, 
    const Eigen::MatrixBase<TangentVectorType>& v) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  pinocchio::forwardKinematics(model_, data_, q, v);
  pinocchio::updateFramePlacements(model_, data_);
  pinocchio::computeForwardKinematicsDerivatives(model_, data_, q, v, 
                                                 Eigen::VectorXd::Zero(dimv_));
}


template <typename VectorType>
inline void Robot::computeBaumgarteResidual(
    const ContactStatus& contact_status, const double time_step,
    const Eigen::MatrixBase<VectorType>& baumgarte_residual) const {
  assert(time_step > 0);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_[i].computeBaumgarteResidual(
          model_, data_, time_step,
          (const_cast<Eigen::MatrixBase<VectorType>&>(baumgarte_residual))
              .template segment<3>(3*num_active_contacts));
      ++num_active_contacts;
    }
  }
}


template <typename VectorType>
inline void Robot::computeBaumgarteResidual(
    const ContactStatus& contact_status, const double coeff, 
    const double time_step,
    const Eigen::MatrixBase<VectorType>& baumgarte_residual) const {
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_[i].computeBaumgarteResidual(
          model_, data_, coeff, time_step,
          (const_cast<Eigen::MatrixBase<VectorType>&>(baumgarte_residual))
              .template segment<3>(3*num_active_contacts));
      ++num_active_contacts;
    }
  }
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void Robot::computeBaumgarteDerivatives(
    const ContactStatus& contact_status, const double time_step,
    const Eigen::MatrixBase<MatrixType1>& baumgarte_partial_dq, 
    const Eigen::MatrixBase<MatrixType2>& baumgarte_partial_dv, 
    const Eigen::MatrixBase<MatrixType3>& baumgarte_partial_da) {
  assert(baumgarte_partial_dq.cols() == dimv_);
  assert(baumgarte_partial_dv.cols() == dimv_);
  assert(baumgarte_partial_da.cols() == dimv_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_[i].computeBaumgarteDerivatives(
          model_, data_, time_step,
          (const_cast<Eigen::MatrixBase<MatrixType1>&>(baumgarte_partial_dq))
              .block(3*num_active_contacts, 0, 3, dimv_),
          (const_cast<Eigen::MatrixBase<MatrixType2>&>(baumgarte_partial_dv))
              .block(3*num_active_contacts, 0, 3, dimv_),
          (const_cast<Eigen::MatrixBase<MatrixType3>&>(baumgarte_partial_da))
              .block(3*num_active_contacts, 0, 3, dimv_));
      ++num_active_contacts;
    }
  }
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void Robot::computeBaumgarteDerivatives(
    const ContactStatus& contact_status, const double coeff, 
    const double time_step,
    const Eigen::MatrixBase<MatrixType1>& baumgarte_partial_dq, 
    const Eigen::MatrixBase<MatrixType2>& baumgarte_partial_dv, 
    const Eigen::MatrixBase<MatrixType3>& baumgarte_partial_da) {
  assert(baumgarte_partial_dq.cols() == dimv_);
  assert(baumgarte_partial_dv.cols() == dimv_);
  assert(baumgarte_partial_da.cols() == dimv_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_[i].computeBaumgarteDerivatives(
          model_, data_, coeff, time_step,
          (const_cast<Eigen::MatrixBase<MatrixType1>&>(baumgarte_partial_dq))
              .block(3*num_active_contacts, 0, 3, dimv_),
          (const_cast<Eigen::MatrixBase<MatrixType2>&>(baumgarte_partial_dv))
              .block(3*num_active_contacts, 0, 3, dimv_),
          (const_cast<Eigen::MatrixBase<MatrixType3>&>(baumgarte_partial_da))
              .block(3*num_active_contacts, 0, 3, dimv_));
      ++num_active_contacts;
    }
  }
}


template <typename VectorType>
inline void Robot::computeContactVelocityResidual(
    const ContactStatus& contact_status, 
    const Eigen::MatrixBase<VectorType>& velocity_residual) const {
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_[i].computeContactVelocityResidual(
          model_, data_,
          (const_cast<Eigen::MatrixBase<VectorType>&>(velocity_residual))
              .template segment<3>(3*num_active_contacts));
      ++num_active_contacts;
    }
  }
}


template <typename MatrixType1, typename MatrixType2>
inline void Robot::computeContactVelocityDerivatives(
    const ContactStatus& contact_status, 
    const Eigen::MatrixBase<MatrixType1>& velocity_partial_dq, 
    const Eigen::MatrixBase<MatrixType2>& velocity_partial_dv) {
  assert(velocity_partial_dq.cols() == dimv_);
  assert(velocity_partial_dv.cols() == dimv_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_[i].computeContactVelocityDerivatives(
          model_, data_, 
          (const_cast<Eigen::MatrixBase<MatrixType1>&>(velocity_partial_dq))
              .block(3*num_active_contacts, 0, 3, dimv_),
          (const_cast<Eigen::MatrixBase<MatrixType2>&>(velocity_partial_dv))
              .block(3*num_active_contacts, 0, 3, dimv_));
      ++num_active_contacts;
    }
  }
}


template <typename VectorType>
inline void Robot::computeContactResidual(
    const ContactStatus& contact_status, 
    const Eigen::MatrixBase<VectorType>& contact_residual) const {
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_[i].computeContactResidual(
          model_, data_,
          (const_cast<Eigen::MatrixBase<VectorType>&>(contact_residual))
              .template segment<3>(3*num_active_contacts));
      ++num_active_contacts;
    }
  }
}


template <typename VectorType>
inline void Robot::computeContactResidual(
    const ContactStatus& contact_status, const double coeff,
    const Eigen::MatrixBase<VectorType>& contact_residual) const {
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_[i].computeContactResidual(
          model_, data_, coeff,
          (const_cast<Eigen::MatrixBase<VectorType>&>(contact_residual))
              .template segment<3>(3*num_active_contacts));
      ++num_active_contacts;
    }
  }
}


template <typename MatrixType>
inline void Robot::computeContactDerivative(
    const ContactStatus& contact_status,
    const Eigen::MatrixBase<MatrixType>& contact_partial_dq) {
  assert(contact_partial_dq.cols() == dimv_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_[i].computeContactDerivative(
          model_, data_, 
          (const_cast<Eigen::MatrixBase<MatrixType>&>(contact_partial_dq))
              .block(3*num_active_contacts, 0, 3, dimv_));
      ++num_active_contacts;
    }
  }
}


template <typename MatrixType>
inline void Robot::computeContactDerivative(
    const ContactStatus& contact_status, const double coeff,
    const Eigen::MatrixBase<MatrixType>& contact_partial_dq) {
  assert(contact_partial_dq.cols() == dimv_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_[i].computeContactDerivative(
          model_, data_, coeff,
          (const_cast<Eigen::MatrixBase<MatrixType>&>(contact_partial_dq))
              .block(3*num_active_contacts, 0, 3, dimv_));
      ++num_active_contacts;
    }
  }
}


inline void Robot::setContactPoints(
    const std::vector<Eigen::Vector3d>& contact_points) {
  assert(contact_points.size() == point_contacts_.size());
  for (int i=0; i<point_contacts_.size(); ++i) {
    point_contacts_[i].setContactPoint(contact_points[i]);
  }
}


inline void Robot::setContactPointsByCurrentKinematics() {
  for (int i=0; i<point_contacts_.size(); ++i) {
    point_contacts_[i].setContactPointByCurrentKinematics(data_);
  }
}


inline void Robot::setFrictionCoefficient(
    const std::vector<double>& friction_coefficient) {
  for (int i=0; i<point_contacts_.size(); ++i) {
    point_contacts_[i].setFrictionCoefficient(friction_coefficient[i]);
  }
}


inline double Robot::frictionCoefficient(const int contact_index) const {
  assert(contact_index >= 0);
  assert(contact_index < point_contacts_.size());
  try {
    if (point_contacts_.empty()) {
      throw std::runtime_error(
          "invalid function call: robot has no point contacts!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  return point_contacts_[contact_index].frictionCoefficient();
}


inline void Robot::setRestitutionCoefficient(
    const std::vector<double>& restitution_coefficient) {
  for (int i=0; i<point_contacts_.size(); ++i) {
    point_contacts_[i].setRestitutionCoefficient(restitution_coefficient[i]);
  }
}


inline double Robot::restitutionCoefficient(const int contact_index) const {
  assert(contact_index >= 0);
  assert(contact_index < point_contacts_.size());
  try {
    if (point_contacts_.empty()) {
      throw std::runtime_error(
          "invalid function call: robot has no point contacts!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  return point_contacts_[contact_index].restitutionCoefficient();
}


inline void Robot::setContactForces(const ContactStatus& contact_status, 
                                    const std::vector<Eigen::Vector3d>& f) {
  assert(f.size() == point_contacts_.size());
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_[i].computeJointForceFromContactForce(f[i], fjoint_);
      ++num_active_contacts;
    }
    else {
      point_contacts_[i].computeJointForceFromContactForce(
          Eigen::Vector3d::Zero(), fjoint_);
    }
  }
}


template <typename ConfigVectorType, typename TangentVectorType1, 
          typename TangentVectorType2, typename TangentVectorType3>
inline void Robot::RNEA(const Eigen::MatrixBase<ConfigVectorType>& q, 
                        const Eigen::MatrixBase<TangentVectorType1>& v, 
                        const Eigen::MatrixBase<TangentVectorType2>& a, 
                        const Eigen::MatrixBase<TangentVectorType3>& tau) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(tau.size() == dimv_);
  if (point_contacts_.empty()) {
    const_cast<Eigen::MatrixBase<TangentVectorType3>&>(tau)
        = pinocchio::rnea(model_, data_, q, v, a);
  }
  else {
    const_cast<Eigen::MatrixBase<TangentVectorType3>&>(tau)
        = pinocchio::rnea(model_, data_, q, v, a, fjoint_);
  }
}


template <typename ConfigVectorType, typename TangentVectorType1, 
          typename TangentVectorType2, typename MatrixType1, 
          typename MatrixType2, typename MatrixType3>
inline void Robot::RNEADerivatives(
    const Eigen::MatrixBase<ConfigVectorType>& q, 
    const Eigen::MatrixBase<TangentVectorType1>& v, 
    const Eigen::MatrixBase<TangentVectorType2>& a, 
    const Eigen::MatrixBase<MatrixType1>& dRNEA_partial_dq, 
    const Eigen::MatrixBase<MatrixType2>& dRNEA_partial_dv, 
    const Eigen::MatrixBase<MatrixType3>& dRNEA_partial_da) {
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
      pinocchio::computeRNEADerivatives(
          model_, data_, q, v, a, 
          const_cast<Eigen::MatrixBase<MatrixType1>&>(dRNEA_partial_dq),
          const_cast<Eigen::MatrixBase<MatrixType2>&>(dRNEA_partial_dv),
          const_cast<Eigen::MatrixBase<MatrixType3>&>(dRNEA_partial_da));
  }
  else {
      pinocchio::computeRNEADerivatives(
          model_, data_, q, v, a, fjoint_,
          const_cast<Eigen::MatrixBase<MatrixType1>&>(dRNEA_partial_dq),
          const_cast<Eigen::MatrixBase<MatrixType2>&>(dRNEA_partial_dv),
          const_cast<Eigen::MatrixBase<MatrixType3>&>(dRNEA_partial_da));
  }
  (const_cast<Eigen::MatrixBase<MatrixType3>&>(dRNEA_partial_da)) 
      .template triangularView<Eigen::StrictlyLower>() 
      = (const_cast<Eigen::MatrixBase<MatrixType3>&>(dRNEA_partial_da)).transpose()
          .template triangularView<Eigen::StrictlyLower>();
}


template <typename ConfigVectorType, typename TangentVectorType1, 
          typename TangentVectorType2>
inline void Robot::RNEAImpulse(
    const Eigen::MatrixBase<ConfigVectorType>& q, 
    const Eigen::MatrixBase<TangentVectorType1>& dv,
    const Eigen::MatrixBase<TangentVectorType2>& res) {
  assert(q.size() == dimq_);
  assert(dv.size() == dimv_);
  assert(res.size() == dimv_);
  const_cast<Eigen::MatrixBase<TangentVectorType2>&>(res)
      = pinocchio::rnea(impulse_model_, impulse_data_, q, 
                        Eigen::VectorXd::Zero(dimv_),  dv, fjoint_);
}


template <typename ConfigVectorType, typename TangentVectorType, 
          typename MatrixType1, typename MatrixType2>
inline void Robot::RNEAImpulseDerivatives(
    const Eigen::MatrixBase<ConfigVectorType>& q, 
    const Eigen::MatrixBase<TangentVectorType>& dv, 
    const Eigen::MatrixBase<MatrixType1>& dRNEA_partial_dq, 
    const Eigen::MatrixBase<MatrixType2>& dRNEA_partial_ddv) {
  assert(q.size() == dimq_);
  assert(dv.size() == dimv_);
  assert(dRNEA_partial_dq.cols() == dimv_);
  assert(dRNEA_partial_dq.rows() == dimv_);
  assert(dRNEA_partial_ddv.cols() == dimv_);
  assert(dRNEA_partial_ddv.rows() == dimv_);
  pinocchio::computeRNEADerivatives(
      impulse_model_, impulse_data_, q, Eigen::VectorXd::Zero(dimv_), dv, 
      fjoint_, const_cast<Eigen::MatrixBase<MatrixType1>&>(dRNEA_partial_dq),
      dimpulse_dv_,
      const_cast<Eigen::MatrixBase<MatrixType2>&>(dRNEA_partial_ddv));
  (const_cast<Eigen::MatrixBase<MatrixType2>&>(dRNEA_partial_ddv)) 
      .template triangularView<Eigen::StrictlyLower>() 
      = (const_cast<Eigen::MatrixBase<MatrixType2>&>(dRNEA_partial_ddv)).transpose()
          .template triangularView<Eigen::StrictlyLower>();
}


template <typename MatrixType>
inline void Robot::dRNEAPartialdFext(
    const ContactStatus& contact_status,
    const Eigen::MatrixBase<MatrixType>& dRNEA_partial_dfext) {
  assert(dRNEA_partial_dfext.rows() == dimv_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_[i].getContactJacobian(
          model_, data_,  -1, 
          (const_cast<Eigen::MatrixBase<MatrixType>&>(dRNEA_partial_dfext))
              .block(0, 3*num_active_contacts, dimv_, 3), true);
      ++num_active_contacts;
    }
  }
}


template <typename ConfigVectorType, typename TangentVectorType1, 
          typename TangentVectorType2, typename TangentVectorType3,
          typename TangentVectorType4>
inline void Robot::stateEquation(
    const Eigen::MatrixBase<ConfigVectorType>& q, 
    const Eigen::MatrixBase<TangentVectorType1>& v, 
    const Eigen::MatrixBase<TangentVectorType2>& tau, 
    const Eigen::MatrixBase<TangentVectorType3>& dq, 
    const Eigen::MatrixBase<TangentVectorType4>& dv) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(tau.size() == dimv_);
  assert(dq.size() == dimv_);
  assert(dv.size() == dimv_);
  const_cast<Eigen::MatrixBase<TangentVectorType3>&> (dq) = v;
  if (point_contacts_.empty()) {
    const_cast<Eigen::MatrixBase<TangentVectorType3>&> (dv)
        = pinocchio::aba(model_, data_, q, v, tau);
  }
  else {
    const_cast<Eigen::MatrixBase<TangentVectorType3>&> (dv)
        = pinocchio::aba(model_, data_, q, v, tau, fjoint_);
  }
}


template <typename ConfigVectorType>
inline void Robot::generateFeasibleConfiguration(
    const Eigen::MatrixBase<ConfigVectorType>& q) const {
  assert(q.size() == dimq_);
  Eigen::VectorXd q_min = model_.lowerPositionLimit;
  Eigen::VectorXd q_max = model_.upperPositionLimit;
  if (floating_base_.has_floating_base()) {
    q_min.head(7) = - Eigen::VectorXd::Ones(7);
    q_max.head(7) = Eigen::VectorXd::Ones(7);
  }
  const_cast<Eigen::MatrixBase<ConfigVectorType>&> (q) 
      = pinocchio::randomConfiguration(model_, q_min, q_max);
}


template <typename ConfigVectorType>
inline void Robot::normalizeConfiguration(
    const Eigen::MatrixBase<ConfigVectorType>& q) const {
  assert(q.size() == dimq_);
  if (floating_base_.has_floating_base()) {
    if (q.template segment<4>(3).squaredNorm() 
          <= std::numeric_limits<double>::epsilon()) {
      (const_cast<Eigen::MatrixBase<ConfigVectorType>&> (q)).coeffRef(3) = 1;
    }
    pinocchio::normalize(model_, 
                         const_cast<Eigen::MatrixBase<ConfigVectorType>&>(q));
  }
}


inline Eigen::VectorXd Robot::jointEffortLimit() const {
  return joint_effort_limit_;
}


inline Eigen::VectorXd Robot::jointVelocityLimit() const {
  return joint_velocity_limit_;
}


inline Eigen::VectorXd Robot::lowerJointPositionLimit() const {
  return lower_joint_position_limit_;
}


inline Eigen::VectorXd Robot::upperJointPositionLimit() const {
  return upper_joint_position_limit_;
}


inline int Robot::dimq() const {
  return dimq_;
}


inline int Robot::dimv() const {
  return dimv_;
}


inline int Robot::max_dimf() const {
  return max_dimf_;
}


inline int Robot::dim_passive() const {
  return floating_base_.dim_passive();
}


inline bool Robot::has_floating_base() const {
  return floating_base_.has_floating_base();
}


inline int Robot::max_point_contacts() const {
  return point_contacts_.size();
}


inline void Robot::initializeJointLimits() {
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

#endif // IDOCP_ROBOT_HXX_ 