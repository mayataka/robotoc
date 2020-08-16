#include "idocp/robot/point_contact.hpp"

#include "assert.h"


namespace idocp {

template <typename MatrixType>
void PointContact::getContactJacobian(
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
void PointContact::getContactJacobian(
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
void PointContact::computeBaumgarteResidual(
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
void PointContact::computeBaumgarteResidual(
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
void PointContact::computeBaumgarteDerivatives(
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
void PointContact::computeBaumgarteDerivatives(
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

} // namespace idocp