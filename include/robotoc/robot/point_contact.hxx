#ifndef ROBOTOC_POINT_CONTACT_HXX_
#define ROBOTOC_POINT_CONTACT_HXX_

#include "robotoc/robot/point_contact.hpp"

#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/kinematics-derivatives.hpp"
#include "pinocchio/algorithm/frames-derivatives.hpp"

#include <cassert>


namespace robotoc {

template <typename VectorType1, typename VectorType2>
inline void PointContact::computeBaumgarteResidual(
    const pinocchio::Model& model, const pinocchio::Data& data, 
    const Eigen::MatrixBase<VectorType1>& desired_contact_position,
    const Eigen::MatrixBase<VectorType2>& baumgarte_residual) const {
  assert(desired_contact_position.size() == 3);
  assert(baumgarte_residual.size() == 3);
  const_cast<Eigen::MatrixBase<VectorType2>&> (baumgarte_residual).noalias()
      = pinocchio::getFrameClassicalAcceleration(model, data, contact_frame_id_, 
                                                 pinocchio::LOCAL).linear();
  (const_cast<Eigen::MatrixBase<VectorType2>&> (baumgarte_residual)).noalias()
      += info_.baumgarte_velocity_gain 
          * pinocchio::getFrameVelocity(model, data, contact_frame_id_, 
                                              pinocchio::LOCAL).linear();
  (const_cast<Eigen::MatrixBase<VectorType2>&> (baumgarte_residual)).noalias()
      += info_.baumgarte_position_gain
          * (data.oMf[contact_frame_id_].translation()-desired_contact_position);
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void PointContact::computeBaumgarteDerivatives(
    const pinocchio::Model& model, pinocchio::Data& data, 
    const Eigen::MatrixBase<MatrixType1>& baumgarte_partial_dq, 
    const Eigen::MatrixBase<MatrixType2>& baumgarte_partial_dv, 
    const Eigen::MatrixBase<MatrixType3>& baumgarte_partial_da) {
  assert(baumgarte_partial_dq.cols() == dimv_);
  assert(baumgarte_partial_dv.cols() == dimv_);
  assert(baumgarte_partial_da.cols() == dimv_);
  assert(baumgarte_partial_dq.rows() == 3);
  assert(baumgarte_partial_dv.rows() == 3);
  assert(baumgarte_partial_da.rows() == 3);
  pinocchio::getFrameAccelerationDerivatives(model, data, contact_frame_id_, 
                                             pinocchio::LOCAL,
                                             frame_v_partial_dq_, 
                                             frame_a_partial_dq_, 
                                             frame_a_partial_dv_, 
                                             frame_a_partial_da_);
  // Skew matrices and LOCAL frame Jacobian are needed to convert the 
  // frame acceleration derivatives into the "classical" acceleration 
  // derivatives.
  pinocchio::getFrameJacobian(model, data, contact_frame_id_,  
                              pinocchio::LOCAL, J_frame_);
  v_frame_ = pinocchio::getFrameVelocity(model, data, contact_frame_id_, 
                                         pinocchio::LOCAL);
  pinocchio::skew(v_frame_.linear(), v_linear_skew_);
  pinocchio::skew(v_frame_.angular(), v_angular_skew_);
  const_cast<Eigen::MatrixBase<MatrixType1>&> (baumgarte_partial_dq)
      = frame_a_partial_dq_.template topRows<3>();
  const_cast<Eigen::MatrixBase<MatrixType1>&> (baumgarte_partial_dq).noalias()
      += v_angular_skew_ * frame_v_partial_dq_.template topRows<3>();
  const_cast<Eigen::MatrixBase<MatrixType1>&> (baumgarte_partial_dq).noalias()
      -= v_linear_skew_ * frame_v_partial_dq_.template bottomRows<3>();
  const_cast<Eigen::MatrixBase<MatrixType2>&> (baumgarte_partial_dv)
      = frame_a_partial_dv_.template topRows<3>();
  const_cast<Eigen::MatrixBase<MatrixType2>&> (baumgarte_partial_dv).noalias()
      += v_angular_skew_ * J_frame_.template topRows<3>();
  const_cast<Eigen::MatrixBase<MatrixType2>&> (baumgarte_partial_dv).noalias()
      -= v_linear_skew_ * J_frame_.template bottomRows<3>();
  const_cast<Eigen::MatrixBase<MatrixType3>&> (baumgarte_partial_da)
      = frame_a_partial_da_.template topRows<3>();
  (const_cast<Eigen::MatrixBase<MatrixType1>&> (baumgarte_partial_dq)).noalias()
      += info_.baumgarte_velocity_gain
          * frame_v_partial_dq_.template topRows<3>();
  (const_cast<Eigen::MatrixBase<MatrixType2>&> (baumgarte_partial_dv)).noalias() 
      += info_.baumgarte_velocity_gain 
          * frame_a_partial_da_.template topRows<3>();
  (const_cast<Eigen::MatrixBase<MatrixType1>&> (baumgarte_partial_dq)).noalias()
      += info_.baumgarte_position_gain * data.oMf[contact_frame_id_].rotation()
                                       * J_frame_.template topRows<3>();
}


template <typename VectorType>
inline void PointContact::computeContactVelocityResidual(
    const pinocchio::Model& model, const pinocchio::Data& data, 
    const Eigen::MatrixBase<VectorType>& velocity_residual) const {
  assert(velocity_residual.size() == 3);
  const_cast<Eigen::MatrixBase<VectorType>&> (velocity_residual).noalias()
      = pinocchio::getFrameVelocity(model, data, contact_frame_id_, 
                                    pinocchio::LOCAL).linear();
}


template <typename MatrixType1, typename MatrixType2>
inline void PointContact::computeContactVelocityDerivatives(
    const pinocchio::Model& model, pinocchio::Data& data,
    const Eigen::MatrixBase<MatrixType1>& velocity_partial_dq, 
    const Eigen::MatrixBase<MatrixType2>& velocity_partial_dv) {
  assert(velocity_partial_dq.cols() == dimv_);
  assert(velocity_partial_dv.cols() == dimv_);
  assert(velocity_partial_dq.rows() == 3);
  assert(velocity_partial_dv.rows() == 3);
  pinocchio::getFrameVelocityDerivatives(model, data, contact_frame_id_, 
                                         pinocchio::LOCAL,
                                         frame_v_partial_dq_, 
                                         frame_a_partial_da_);
  (const_cast<Eigen::MatrixBase<MatrixType1>&> (velocity_partial_dq)).noalias()
      = frame_v_partial_dq_.template topRows<3>();
  (const_cast<Eigen::MatrixBase<MatrixType2>&> (velocity_partial_dv)).noalias() 
      = frame_a_partial_da_.template topRows<3>();
}


template <typename VectorType1, typename VectorType2>
inline void PointContact::computeContactPositionResidual(
    const pinocchio::Model& model, const pinocchio::Data& data, 
    const Eigen::MatrixBase<VectorType1>& desired_contact_position,
    const Eigen::MatrixBase<VectorType2>& position_residual) const {
  assert(desired_contact_position.size() == 3);
  assert(position_residual.size() == 3);
  (const_cast<Eigen::MatrixBase<VectorType2>&> (position_residual))
      = (data.oMf[contact_frame_id_].translation()-desired_contact_position);
}


template <typename MatrixType>
inline void PointContact::computeContactPositionDerivative(
    const pinocchio::Model& model, pinocchio::Data& data, 
    const Eigen::MatrixBase<MatrixType>& position_partial_dq) {
  assert(position_partial_dq.cols() == dimv_);
  assert(position_partial_dq.rows() == 3);
  pinocchio::getFrameJacobian(model, data, contact_frame_id_,  
                              pinocchio::LOCAL, J_frame_);
  (const_cast<Eigen::MatrixBase<MatrixType>&> (position_partial_dq)).noalias()
      = data.oMf[contact_frame_id_].rotation() * J_frame_.template topRows<3>();
}

} // namespace robotoc

#endif // ROBOTOC_POINT_CONTACT_HXX_