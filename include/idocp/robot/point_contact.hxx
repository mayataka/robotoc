#ifndef IDOCP_POINT_CONTACT_HXX_
#define IDOCP_POINT_CONTACT_HXX_

#include "assert.h"


namespace idocp {

inline void PointContact::resetContactPoint(
    const Eigen::Vector3d& contact_point) {
  contact_point_ = contact_point;
}


inline void PointContact::resetContactPointByCurrentKinematics(
    const pinocchio::Data& data) {
  contact_point_ = data.oMf[contact_frame_id_].translation();
}


inline void PointContact::computeJointForceFromContactForce(
    const Eigen::Vector3d& contact_force, 
    pinocchio::container::aligned_vector<pinocchio::Force>& joint_forces) const {
  joint_forces[parent_joint_id_] 
      = jXf_.act(pinocchio::Force(contact_force, Eigen::Vector3d::Zero()));
}


template <typename MatrixType>
inline void PointContact::getContactJacobian(
    const pinocchio::Model& model, pinocchio::Data& data, 
    const Eigen::MatrixBase<MatrixType>& Jacobian, const bool transpose) {
  pinocchio::getFrameJacobian(model, data, contact_frame_id_,  
                              pinocchio::LOCAL, J_frame_);
  if (transpose) {
      assert(Jacobian.rows() == dimv_);
      assert(Jacobian.cols() == 3);
      const_cast<Eigen::MatrixBase<MatrixType>&> (Jacobian)
          = J_frame_.topRows<3>().transpose(); 
  } 
  else {
      assert(Jacobian.rows() == 3);
      assert(Jacobian.cols() == dimv_);
      const_cast<Eigen::MatrixBase<MatrixType>&> (Jacobian)
          = J_frame_.topRows<3>(); 
  }
}


template <typename MatrixType>
inline void PointContact::getContactJacobian(
    const pinocchio::Model& model, pinocchio::Data& data, 
    const double coeff, const Eigen::MatrixBase<MatrixType>& Jacobian, 
    const bool transpose) {
  pinocchio::getFrameJacobian(model, data, contact_frame_id_,  
                              pinocchio::LOCAL, J_frame_);
  pinocchio::getFrameJacobian(model, data, contact_frame_id_,  
                              pinocchio::LOCAL, J_frame_);
  if (transpose) {
      assert(Jacobian.rows() == dimv_);
      assert(Jacobian.cols() == 3);
      const_cast<Eigen::MatrixBase<MatrixType>&> (Jacobian)
          = coeff * J_frame_.topRows<3>().transpose(); 
  } 
  else {
      assert(Jacobian.rows() == 3);
      assert(Jacobian.cols() == dimv_);
      const_cast<Eigen::MatrixBase<MatrixType>&> (Jacobian)
          = coeff * J_frame_.topRows<3>(); 
  }
}


template <typename VectorType>
inline void PointContact::computeBaumgarteResidual(
    const pinocchio::Model& model, const pinocchio::Data& data, 
    const Eigen::MatrixBase<VectorType>& baumgarte_residual) const {
  assert(baumgarte_residual.size() == 3);
  const_cast<Eigen::MatrixBase<VectorType>&> (baumgarte_residual)
      = pinocchio::getFrameClassicalAcceleration(model, data, 
                                                  contact_frame_id_, 
                                                  pinocchio::LOCAL).linear();
  if (baumgarte_weight_on_velocity_ != 0.) {
    (const_cast<Eigen::MatrixBase<VectorType>&> (baumgarte_residual)).noalias()
        += baumgarte_weight_on_velocity_ 
              * pinocchio::getFrameVelocity(model, data, contact_frame_id_, 
                                            pinocchio::LOCAL).linear();
  }
  if (baumgarte_weight_on_position_ != 0.) {
    (const_cast<Eigen::MatrixBase<VectorType>&> (baumgarte_residual)).noalias()
        += baumgarte_weight_on_position_
              * (data.oMf[contact_frame_id_].translation()-contact_point_);
  }
}


template <typename VectorType>
inline void PointContact::computeBaumgarteResidual(
    const pinocchio::Model& model, const pinocchio::Data& data, 
    const double coeff, 
    const Eigen::MatrixBase<VectorType>& baumgarte_residual) const {
  assert(baumgarte_residual.size() == 3);
  const_cast<Eigen::MatrixBase<VectorType>&> (baumgarte_residual)
      = coeff * pinocchio::getFrameClassicalAcceleration(
                    model, data, contact_frame_id_, pinocchio::LOCAL).linear();
  if (baumgarte_weight_on_velocity_ != 0.) {
    (const_cast<Eigen::MatrixBase<VectorType>&> (baumgarte_residual)).noalias()
        += coeff * baumgarte_weight_on_velocity_ 
                 * pinocchio::getFrameVelocity(model, data, contact_frame_id_, 
                                               pinocchio::LOCAL).linear();
  }
  if (baumgarte_weight_on_position_ != 0.) {
    (const_cast<Eigen::MatrixBase<VectorType>&> (baumgarte_residual)).noalias()
        += coeff * baumgarte_weight_on_position_
                 * (data.oMf[contact_frame_id_].translation()-contact_point_);
  }
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void PointContact::computeBaumgarteDerivatives(
    const pinocchio::Model& model, pinocchio::Data& data, 
    const Eigen::MatrixBase<MatrixType1>& baumgarte_partial_dq, 
    const Eigen::MatrixBase<MatrixType2>& baumgarte_partial_dv, 
    const Eigen::MatrixBase<MatrixType3>& baumgarte_partial_da) {
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
      = frame_a_partial_dq_.template topRows<3>()
          + v_angular_skew_ * frame_v_partial_dq_.template topRows<3>()
          + v_linear_skew_ * frame_v_partial_dq_.template bottomRows<3>();
  const_cast<Eigen::MatrixBase<MatrixType2>&> (baumgarte_partial_dv)
      = frame_a_partial_dv_.template topRows<3>()
          + v_angular_skew_ * J_frame_.template topRows<3>()
          + v_linear_skew_ * J_frame_.template bottomRows<3>();
  const_cast<Eigen::MatrixBase<MatrixType3>&> (baumgarte_partial_da)
      = frame_a_partial_da_.template topRows<3>();
  if (baumgarte_weight_on_velocity_ != 0.) {
    (const_cast<Eigen::MatrixBase<MatrixType1>&> (baumgarte_partial_dq)).noalias()
        += baumgarte_weight_on_velocity_ 
            * frame_v_partial_dq_.template topRows<3>();
    (const_cast<Eigen::MatrixBase<MatrixType2>&> (baumgarte_partial_dv)).noalias() 
        += baumgarte_weight_on_velocity_ 
            * frame_a_partial_da_.template topRows<3>();
  }
  if (baumgarte_weight_on_position_ != 0.) {
    (const_cast<Eigen::MatrixBase<MatrixType1>&> (baumgarte_partial_dq)).noalias()
        += baumgarte_weight_on_position_ 
            * data.oMf[contact_frame_id_].rotation()
            * J_frame_.template topRows<3>();
  }
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void PointContact::computeBaumgarteDerivatives(
    const pinocchio::Model& model, pinocchio::Data& data, const double coeff,
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
      = coeff * frame_a_partial_dq_.template topRows<3>()
          + coeff * v_angular_skew_ * frame_v_partial_dq_.template topRows<3>()
          + coeff * v_linear_skew_ * frame_v_partial_dq_.template bottomRows<3>();
  const_cast<Eigen::MatrixBase<MatrixType2>&> (baumgarte_partial_dv) 
      = coeff * frame_a_partial_dv_.template topRows<3>()
          + coeff * v_angular_skew_ * J_frame_.template topRows<3>()
          + coeff * v_linear_skew_ * J_frame_.template bottomRows<3>();
  const_cast<Eigen::MatrixBase<MatrixType3>&> (baumgarte_partial_da)
      = coeff * frame_a_partial_da_.template topRows<3>();
  if (baumgarte_weight_on_velocity_ != 0.) {
    (const_cast<Eigen::MatrixBase<MatrixType1>&> (baumgarte_partial_dq)).noalias()
        += coeff * baumgarte_weight_on_velocity_ 
            * frame_v_partial_dq_.template topRows<3>();
    (const_cast<Eigen::MatrixBase<MatrixType2>&> (baumgarte_partial_dv)).noalias() 
        += coeff * baumgarte_weight_on_velocity_ 
            * frame_a_partial_da_.template topRows<3>();
  }
  if (baumgarte_weight_on_position_ != 0.) {
    (const_cast<Eigen::MatrixBase<MatrixType1>&> (baumgarte_partial_dq)).noalias()
        += coeff * baumgarte_weight_on_position_ 
            * data.oMf[contact_frame_id_].rotation()
            * J_frame_.template topRows<3>();
  }
}


inline void PointContact::activate() {
  is_active_ = true;
}


inline void PointContact::deactivate() {
  is_active_ = false;
}


inline bool PointContact::isActive() const {
  return is_active_;
}


inline int PointContact::contact_frame_id() const {
  return contact_frame_id_;
}


inline int PointContact::parent_joint_id() const {
  return parent_joint_id_;
}


inline double PointContact::baumgarte_weight_on_velocity() const {
  return baumgarte_weight_on_velocity_;
}


inline double PointContact::baumgarte_weight_on_position() const {
  return baumgarte_weight_on_position_;
}


inline Eigen::Vector3d PointContact::contact_point() const {
  return contact_point_;
}

} // namespace idocp

#endif // IDOCP_POINT_CONTACT_HXX_