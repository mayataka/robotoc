#ifndef IDOCP_POINT_CONTACT_HXX_
#define IDOCP_POINT_CONTACT_HXX_

#include "idocp/robot/point_contact.hpp"

#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/kinematics-derivatives.hpp"
#include "pinocchio/algorithm/frames-derivatives.hpp"

#include <stdexcept>
#include <cassert>


namespace idocp {

inline PointContact::PointContact(const pinocchio::Model& model, 
                                  const int contact_frame_id,
                                  const double baumgarte_weight_on_velocity,
                                  const double baumgarte_weight_on_position)
  : contact_frame_id_(contact_frame_id),
    parent_joint_id_(model.frames[contact_frame_id_].parent), 
    dimv_(model.nv),
    baumgarte_weight_on_velocity_(baumgarte_weight_on_velocity),
    baumgarte_weight_on_position_(baumgarte_weight_on_position), 
    jXf_(model.frames[contact_frame_id_].placement),
    J_frame_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_v_partial_dq_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_a_partial_dq_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_a_partial_dv_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_a_partial_da_(Eigen::MatrixXd::Zero(6, model.nv)) {
  try {
    if (contact_frame_id_ < 0) {
      throw std::out_of_range(
          "Invalid argument: contact_frame_id must be non-negative!");
    }
    if (baumgarte_weight_on_velocity < 0) {
      throw std::out_of_range(
          "Invalid argument: baumgarte_weight_on_velocity must be non-negative!");
    }
    if (baumgarte_weight_on_position < 0) {
      throw std::out_of_range(
          "Invalid argument: baumgarte_weight_on_position must be non-negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  v_frame_.setZero();
  v_linear_skew_.setZero();
  v_angular_skew_.setZero();
}


inline PointContact::~PointContact() {
}


inline void PointContact::setBaumgarteWeights(
    const double baumgarte_weight_on_velocity,
    const double baumgarte_weight_on_position) {
  try {
    if (baumgarte_weight_on_velocity < 0) {
      throw std::out_of_range(
          "Invalid argument: baumgarte_weight_on_velocity must be non-negative!");
    }
    if (baumgarte_weight_on_position < 0) {
      throw std::out_of_range(
          "Invalid argument: baumgarte_weight_on_position must be non-negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  baumgarte_weight_on_velocity_ = baumgarte_weight_on_velocity;
  baumgarte_weight_on_position_ = baumgarte_weight_on_position;
}


inline void PointContact::computeJointForceFromContactForce(
    const Eigen::Vector3d& contact_force, 
    pinocchio::container::aligned_vector<pinocchio::Force>& joint_forces) const {
  joint_forces[parent_joint_id_] 
      = jXf_.act(pinocchio::Force(contact_force, Eigen::Vector3d::Zero()));
}


template <typename VectorType1, typename VectorType2>
inline void PointContact::computeBaumgarteResidual(
    const pinocchio::Model& model, const pinocchio::Data& data, 
    const Eigen::MatrixBase<VectorType1>& contact_point,
    const Eigen::MatrixBase<VectorType2>& baumgarte_residual) const {
  assert(contact_point.size() == 3);
  assert(baumgarte_residual.size() == 3);
  const_cast<Eigen::MatrixBase<VectorType2>&> (baumgarte_residual).noalias()
      = pinocchio::getFrameClassicalAcceleration(model, data, contact_frame_id_, 
                                                 pinocchio::LOCAL).linear();
  (const_cast<Eigen::MatrixBase<VectorType2>&> (baumgarte_residual)).noalias()
      += baumgarte_weight_on_velocity_ 
          * pinocchio::getFrameVelocity(model, data, contact_frame_id_, 
                                              pinocchio::LOCAL).linear();
  (const_cast<Eigen::MatrixBase<VectorType2>&> (baumgarte_residual)).noalias()
      += baumgarte_weight_on_position_
          * (data.oMf[contact_frame_id_].translation()-contact_point);
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
      += baumgarte_weight_on_velocity_
          * frame_v_partial_dq_.template topRows<3>();
  (const_cast<Eigen::MatrixBase<MatrixType2>&> (baumgarte_partial_dv)).noalias() 
      += baumgarte_weight_on_velocity_ 
          * frame_a_partial_da_.template topRows<3>();
  (const_cast<Eigen::MatrixBase<MatrixType1>&> (baumgarte_partial_dq)).noalias()
      += baumgarte_weight_on_position_ * data.oMf[contact_frame_id_].rotation()
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
    const Eigen::MatrixBase<VectorType1>& contact_point,
    const Eigen::MatrixBase<VectorType2>& contact_residual) const {
  assert(contact_point.size() == 3);
  assert(contact_residual.size() == 3);
  (const_cast<Eigen::MatrixBase<VectorType2>&> (contact_residual))
      = (data.oMf[contact_frame_id_].translation()-contact_point);
}


template <typename MatrixType>
inline void PointContact::computeContactPositionDerivative(
    const pinocchio::Model& model, pinocchio::Data& data, 
    const Eigen::MatrixBase<MatrixType>& contact_partial_dq) {
  assert(contact_partial_dq.cols() == dimv_);
  assert(contact_partial_dq.rows() == 3);
  pinocchio::getFrameJacobian(model, data, contact_frame_id_,  
                              pinocchio::LOCAL, J_frame_);
  (const_cast<Eigen::MatrixBase<MatrixType>&> (contact_partial_dq)).noalias()
      = data.oMf[contact_frame_id_].rotation() * J_frame_.template topRows<3>();
}


inline const Eigen::Vector3d& PointContact::contactPoint(
    const pinocchio::Data& data) const {
  return data.oMf[contact_frame_id_].translation();
}


inline int PointContact::contact_frame_id() const {
  return contact_frame_id_;
}


inline int PointContact::parent_joint_id() const {
  return parent_joint_id_;
}

} // namespace idocp

#endif // IDOCP_POINT_CONTACT_HXX_