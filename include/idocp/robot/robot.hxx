#ifndef IDOCP_ROBOT_HXX_
#define IDOCP_ROBOT_HXX_

#include "idocp/robot/robot.hpp"

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


template <typename VectorType>
inline void Robot::computeBaumgarteResidual(
    const Eigen::MatrixBase<VectorType>& baumgarte_residual) const {
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].computeBaumgarteResidual(
          model_, data_, 
          (const_cast<Eigen::MatrixBase<VectorType>&>(baumgarte_residual))
              .template segment<3>(3*num_active_contacts));
      ++num_active_contacts;
    }
  }
}


template <typename VectorType>
inline void Robot::computeBaumgarteResidual(
    const double coeff, 
    const Eigen::MatrixBase<VectorType>& baumgarte_residual) const {
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].computeBaumgarteResidual(
          model_, data_, coeff, 
          (const_cast<Eigen::MatrixBase<VectorType>&>(baumgarte_residual))
              .template segment<3>(3*num_active_contacts));
      ++num_active_contacts;
    }
  }
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void Robot::computeBaumgarteDerivatives(
    const Eigen::MatrixBase<MatrixType1>& baumgarte_partial_dq, 
    const Eigen::MatrixBase<MatrixType2>& baumgarte_partial_dv, 
    const Eigen::MatrixBase<MatrixType3>& baumgarte_partial_da) {
  assert(baumgarte_partial_dq.cols() == dimv_);
  assert(baumgarte_partial_dv.cols() == dimv_);
  assert(baumgarte_partial_da.cols() == dimv_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].computeBaumgarteDerivatives(
          model_, data_, 
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
    const double coeff, 
    const Eigen::MatrixBase<MatrixType1>& baumgarte_partial_dq, 
    const Eigen::MatrixBase<MatrixType2>& baumgarte_partial_dv, 
    const Eigen::MatrixBase<MatrixType3>& baumgarte_partial_da) {
  assert(baumgarte_partial_dq.cols() == dimv_);
  assert(baumgarte_partial_dv.cols() == dimv_);
  assert(baumgarte_partial_da.cols() == dimv_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].computeBaumgarteDerivatives(
          model_, data_, coeff,
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


inline void Robot::setContactPoints(
    const std::vector<Eigen::Vector3d>& contact_points) {
  for (int i=0; i<point_contacts_.size(); ++i) {
    point_contacts_[i].resetContactPoint(contact_points[i]);
  }
}


inline void Robot::setContactPointsByCurrentKinematics() {
  for (int i=0; i<point_contacts_.size(); ++i) {
    point_contacts_[i].resetContactPointByCurrentKinematics(data_);
  }
}


inline void Robot::setContactStatus(
    const std::vector<bool>& is_each_contact_active) {
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
  if (num_active_contacts_ > 0) {
    has_active_contacts_ = true;
  }
  else {
    has_active_contacts_ = false;
  }
}


template <typename VectorType>
inline void Robot::setContactForces(const Eigen::MatrixBase<VectorType>& f) {
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
      point_contacts_[i].computeJointForceFromContactForce(
          f.template segment<3>(3*num_active_contacts), fjoint_);
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


template <typename MatrixType>
inline void Robot::dRNEAPartialdFext(
    const Eigen::MatrixBase<MatrixType>& dRNEA_partial_dfext) {
  assert(dRNEA_partial_dfext.rows() == dimv_);
  int num_active_contacts = 0;
  for (int i=0; i<point_contacts_.size(); ++i) {
    if (point_contacts_[i].isActive()) {
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


inline int Robot::dimJ() const {
  return dimJ_;
}


inline int Robot::max_dimf() const {
  return max_dimf_;
}


inline int Robot::dimf() const {
  return dimf_;
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


inline bool Robot::has_active_contacts() const {
  return has_active_contacts_;
}


inline int Robot::num_active_point_contacts() const {
  return num_active_contacts_;
}


inline bool Robot::is_contact_active(const int contact_index) const {
  assert(contact_index >= 0);
  assert(contact_index < point_contacts_.size());
  return point_contacts_[contact_index].isActive();
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