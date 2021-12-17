#ifndef ROBOTOC_SURFACE_CONTACT_HXX_
#define ROBOTOC_SURFACE_CONTACT_HXX_

#include "robotoc/robot/surface_contact.hpp"

#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/kinematics-derivatives.hpp"
#include "pinocchio/algorithm/frames-derivatives.hpp"

#include <cassert>


namespace robotoc {

inline void SurfaceContact::computeJointForceFromContactWrench(
    const Vector6d& contact_wrench, 
    pinocchio::container::aligned_vector<pinocchio::Force>& joint_forces) const {
  joint_forces[parent_joint_id_] = jXf_.act(pinocchio::Force(contact_wrench));
}


template <typename VectorType>
inline void SurfaceContact::computeBaumgarteResidual(
    const pinocchio::Model& model, const pinocchio::Data& data, 
    const SE3& cotact_placement, 
    const Eigen::MatrixBase<VectorType>& baumgarte_residual) {
  assert(baumgarte_residual.size() == 6);
  const_cast<Eigen::MatrixBase<VectorType>&> (baumgarte_residual).noalias()
      = pinocchio::getFrameAcceleration(model, data, contact_frame_id_, 
                                        pinocchio::LOCAL).toVector();
  (const_cast<Eigen::MatrixBase<VectorType>&> (baumgarte_residual)).noalias()
      += baumgarte_weight_on_velocity_ 
          * pinocchio::getFrameVelocity(model, data, contact_frame_id_, 
                                        pinocchio::LOCAL).toVector();
  X_diff_ = cotact_placement.inverse() * data.oMf[contact_frame_id_];
  (const_cast<Eigen::MatrixBase<VectorType>&> (baumgarte_residual)).noalias()
      += baumgarte_weight_on_position_ * Log6Map(X_diff_);
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void SurfaceContact::computeBaumgarteDerivatives(
    const pinocchio::Model& model, pinocchio::Data& data, 
    const Eigen::MatrixBase<MatrixType1>& baumgarte_partial_dq, 
    const Eigen::MatrixBase<MatrixType2>& baumgarte_partial_dv, 
    const Eigen::MatrixBase<MatrixType3>& baumgarte_partial_da) {
  assert(baumgarte_partial_dq.cols() == dimv_);
  assert(baumgarte_partial_dv.cols() == dimv_);
  assert(baumgarte_partial_da.cols() == dimv_);
  assert(baumgarte_partial_dq.rows() == 6);
  assert(baumgarte_partial_dv.rows() == 6);
  assert(baumgarte_partial_da.rows() == 6);
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
  const_cast<Eigen::MatrixBase<MatrixType1>&> (baumgarte_partial_dq)
      = frame_a_partial_dq_;
  const_cast<Eigen::MatrixBase<MatrixType1>&> (baumgarte_partial_dv)
      = frame_a_partial_dv_;
  const_cast<Eigen::MatrixBase<MatrixType3>&> (baumgarte_partial_da)
      = frame_a_partial_da_;
  (const_cast<Eigen::MatrixBase<MatrixType1>&> (baumgarte_partial_dq)).noalias()
      += baumgarte_weight_on_velocity_ * frame_v_partial_dq_;
  (const_cast<Eigen::MatrixBase<MatrixType2>&> (baumgarte_partial_dv)).noalias() 
      += baumgarte_weight_on_velocity_ * frame_a_partial_da_;
  computeJLog6Map(X_diff_, Jlog6_);
  (const_cast<Eigen::MatrixBase<MatrixType1>&> (baumgarte_partial_dq)).noalias()
      += baumgarte_weight_on_position_ * Jlog6_ * J_frame_;
}


template <typename VectorType>
inline void SurfaceContact::computeContactVelocityResidual(
    const pinocchio::Model& model, const pinocchio::Data& data, 
    const Eigen::MatrixBase<VectorType>& velocity_residual) const {
  assert(velocity_residual.size() == 6);
  const_cast<Eigen::MatrixBase<VectorType>&> (velocity_residual).noalias()
      = pinocchio::getFrameVelocity(model, data, contact_frame_id_, 
                                    pinocchio::LOCAL).toVector();
}


template <typename MatrixType1, typename MatrixType2>
inline void SurfaceContact::computeContactVelocityDerivatives(
    const pinocchio::Model& model, pinocchio::Data& data,
    const Eigen::MatrixBase<MatrixType1>& velocity_partial_dq, 
    const Eigen::MatrixBase<MatrixType2>& velocity_partial_dv) {
  assert(velocity_partial_dq.cols() == dimv_);
  assert(velocity_partial_dv.cols() == dimv_);
  assert(velocity_partial_dq.rows() == 6);
  assert(velocity_partial_dv.rows() == 6);
  pinocchio::getFrameVelocityDerivatives(model, data, contact_frame_id_, 
                                         pinocchio::LOCAL,
                                         frame_v_partial_dq_, 
                                         frame_a_partial_da_);
  (const_cast<Eigen::MatrixBase<MatrixType1>&> (velocity_partial_dq)).noalias()
      = frame_v_partial_dq_;
  (const_cast<Eigen::MatrixBase<MatrixType2>&> (velocity_partial_dv)).noalias() 
      = frame_a_partial_da_;
}


template <typename VectorType>
inline void SurfaceContact::computeContactPositionResidual(
    const pinocchio::Model& model, const pinocchio::Data& data, 
    const SE3& contact_placement,
    const Eigen::MatrixBase<VectorType>& contact_residual) {
  assert(contact_residual.size() == 6);
  X_diff_ = contact_placement.inverse() * data.oMf[contact_frame_id_];
  (const_cast<Eigen::MatrixBase<VectorType>&> (contact_residual))
      = Log6Map(X_diff_);
}


template <typename MatrixType>
inline void SurfaceContact::computeContactPositionDerivative(
    const pinocchio::Model& model, pinocchio::Data& data, 
    const Eigen::MatrixBase<MatrixType>& contact_partial_dq) {
  assert(contact_partial_dq.cols() == dimv_);
  assert(contact_partial_dq.rows() == 6);
  pinocchio::getFrameJacobian(model, data, contact_frame_id_,  
                              pinocchio::LOCAL, J_frame_);
  computeJLog6Map(X_diff_, Jlog6_);
  (const_cast<Eigen::MatrixBase<MatrixType>&> (contact_partial_dq)).noalias()
      = Jlog6_ * J_frame_;
}


inline const SE3& SurfaceContact::contactPlacement(
    const pinocchio::Data& data) const {
  return data.oMf[contact_frame_id_];
}


inline int SurfaceContact::contact_frame_id() const {
  return contact_frame_id_;
}


inline int SurfaceContact::parent_joint_id() const {
  return parent_joint_id_;
}

} // namespace robotoc

#endif // ROBOTOC_SURFACE_CONTACT_HXX_ 