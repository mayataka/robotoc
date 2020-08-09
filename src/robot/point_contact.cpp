#include "idocp/robot/point_contact.hpp"

#include <assert.h>


namespace idocp {

PointContact::PointContact(const pinocchio::Model& model, 
                           const int contact_frame_id, 
                           const double baumgarte_weight_on_velocity, 
                           const double baumgarte_weight_on_position)
  : is_active_(false),
    contact_frame_id_(contact_frame_id),
    parent_joint_id_(model.frames[contact_frame_id_].parent), 
    dimv_(model.nv),
    baumgarte_weight_on_velocity_(baumgarte_weight_on_velocity),
    baumgarte_weight_on_position_(baumgarte_weight_on_position),
    contact_point_(Eigen::Vector3d::Zero()),
    jXf_(model.frames[contact_frame_id_].placement),
    J_frame_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_v_partial_dq_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_a_partial_dq_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_a_partial_dv_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_a_partial_da_(Eigen::MatrixXd::Zero(6, model.nv)) {
  assert(model.check());
  assert(contact_frame_id_ >= 0);
  assert(baumgarte_weight_on_velocity_ >= 0);
  assert(baumgarte_weight_on_position_ >= 0);
  v_frame_.setZero();
  v_linear_skew_.setZero();
  v_angular_skew_.setZero();
}


PointContact::PointContact() 
  : is_active_(false),
    contact_frame_id_(0),
    parent_joint_id_(0), 
    dimv_(0),
    baumgarte_weight_on_velocity_(0),
    baumgarte_weight_on_position_(0),
    contact_point_(Eigen::Vector3d::Zero()),
    jXf_(),
    J_frame_(),
    frame_v_partial_dq_(),
    frame_a_partial_dq_(),
    frame_a_partial_dv_(),
    frame_a_partial_da_() {
  v_frame_.setZero();
  v_linear_skew_.setZero();
  v_angular_skew_.setZero();
}


PointContact::~PointContact() {
}


void PointContact::resetBaugrarteParameters(
    const double baumgarte_weight_on_velocity, 
    const double baumgarte_weight_on_position) {
  assert(baumgarte_weight_on_velocity >= 0);
  assert(baumgarte_weight_on_position >= 0);
  baumgarte_weight_on_velocity_ = baumgarte_weight_on_velocity;
  baumgarte_weight_on_position_ = baumgarte_weight_on_position;
}


void PointContact::resetContactPoint(const Eigen::Vector3d& contact_point) {
  contact_point_ = contact_point;
}


void PointContact::resetContactPointByCurrentKinematics(
    const pinocchio::Data& data) {
  contact_point_ = data.oMf[contact_frame_id_].translation();
}


void PointContact::computeJointForceFromContactForce(
    const Eigen::Vector3d& contact_force, 
    pinocchio::container::aligned_vector<pinocchio::Force>& joint_forces) const {
  joint_forces[parent_joint_id_] 
      = jXf_.act(pinocchio::Force(contact_force, Eigen::Vector3d::Zero()));
}


void PointContact::getContactJacobian(const pinocchio::Model& model, 
                                      pinocchio::Data& data, 
                                      Eigen::Ref<Matrix3Xd> Jacobian) {
  pinocchio::getFrameJacobian(model, data, contact_frame_id_, pinocchio::LOCAL, 
                              J_frame_);
  assert(Jacobian.cols() == dimv_);
  Jacobian = J_frame_.topRows<3>(); 
}


void PointContact::getContactJacobian(const pinocchio::Model& model, 
                                      pinocchio::Data& data, const double coeff,
                                      Eigen::Ref<Matrix3Xd> Jacobian) {
  pinocchio::getFrameJacobian(model, data, contact_frame_id_, pinocchio::LOCAL, 
                              J_frame_);
  assert(Jacobian.cols() == dimv_);
  Jacobian = coeff * J_frame_.topRows<3>(); 
}


void PointContact::computeBaumgarteResidual(
    const pinocchio::Model& model, const pinocchio::Data& data, 
    Eigen::Ref<Eigen::Vector3d> baumgarte_residual) const {
  assert(baumgarte_residual.size() == 3);
  baumgarte_residual
      = pinocchio::getFrameClassicalAcceleration(model, data, contact_frame_id_, 
                                                 pinocchio::LOCAL).linear();
  if (baumgarte_weight_on_velocity_ != 0.) {
    baumgarte_residual.noalias()
        += baumgarte_weight_on_velocity_ 
              * pinocchio::getFrameVelocity(model, data, contact_frame_id_, 
                                            pinocchio::LOCAL).linear();
  }
  if (baumgarte_weight_on_position_ != 0.) {
    baumgarte_residual.noalias()
        += baumgarte_weight_on_position_
              * (data.oMf[contact_frame_id_].translation()-contact_point_);
  }
}


void PointContact::computeBaumgarteResidual(
    const pinocchio::Model& model, const pinocchio::Data& data, 
    const double coeff, Eigen::Ref<Eigen::Vector3d> baumgarte_residual) const {
  assert(baumgarte_residual.size() == 3);
  baumgarte_residual
      = coeff * pinocchio::getFrameClassicalAcceleration(
                    model, data, contact_frame_id_, pinocchio::LOCAL).linear();
  if (baumgarte_weight_on_velocity_ != 0.) {
    baumgarte_residual.noalias()
        += coeff * baumgarte_weight_on_velocity_ 
                 * pinocchio::getFrameVelocity(model, data, contact_frame_id_, 
                                               pinocchio::LOCAL).linear();
  }
  if (baumgarte_weight_on_position_ != 0.) {
    baumgarte_residual.noalias()
        += coeff * baumgarte_weight_on_position_
                 * (data.oMf[contact_frame_id_].translation()-contact_point_);
  }
}


void PointContact::computeBaumgarteDerivatives(
    const pinocchio::Model& model, pinocchio::Data& data, 
    Eigen::Ref<Matrix3Xd> baumgarte_partial_dq, 
    Eigen::Ref<Matrix3Xd> baumgarte_partial_dv, 
    Eigen::Ref<Matrix3Xd> baumgarte_partial_da) {
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
  baumgarte_partial_dq
      = frame_a_partial_dq_.template topRows<3>()
          + v_angular_skew_ * frame_v_partial_dq_.template topRows<3>()
          + v_linear_skew_ * frame_v_partial_dq_.template bottomRows<3>();
  baumgarte_partial_dv 
      = frame_a_partial_dv_.template topRows<3>()
          + v_angular_skew_ * J_frame_.template topRows<3>()
          + v_linear_skew_ * J_frame_.template bottomRows<3>();
  baumgarte_partial_da 
      = frame_a_partial_da_.template topRows<3>();
  if (baumgarte_weight_on_velocity_ != 0.) {
    baumgarte_partial_dq.noalias()
        += baumgarte_weight_on_velocity_ 
            * frame_v_partial_dq_.template topRows<3>();
    baumgarte_partial_dv.noalias() 
        += baumgarte_weight_on_velocity_ 
            * frame_a_partial_da_.template topRows<3>();
  }
  if (baumgarte_weight_on_position_ != 0.) {
    baumgarte_partial_dq.noalias()
        += baumgarte_weight_on_position_ 
            * data.oMf[contact_frame_id_].rotation()
            * J_frame_.template topRows<3>();
  }
}


void PointContact::computeBaumgarteDerivatives(
    const pinocchio::Model& model, pinocchio::Data& data, const double coeff,
    Eigen::Ref<Matrix3Xd> baumgarte_partial_dq, 
    Eigen::Ref<Matrix3Xd> baumgarte_partial_dv, 
    Eigen::Ref<Matrix3Xd> baumgarte_partial_da) {
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
  baumgarte_partial_dq
      = coeff * frame_a_partial_dq_.template topRows<3>()
          + coeff * v_angular_skew_ * frame_v_partial_dq_.template topRows<3>()
          + coeff * v_linear_skew_ * frame_v_partial_dq_.template bottomRows<3>();
  baumgarte_partial_dv 
      = coeff * frame_a_partial_dv_.template topRows<3>()
          + coeff * v_angular_skew_ * J_frame_.template topRows<3>()
          + coeff * v_linear_skew_ * J_frame_.template bottomRows<3>();
  baumgarte_partial_da 
      = coeff * frame_a_partial_da_.template topRows<3>();
  if (baumgarte_weight_on_velocity_ != 0.) {
    baumgarte_partial_dq.noalias()
        += coeff * baumgarte_weight_on_velocity_ 
            * frame_v_partial_dq_.template topRows<3>();
    baumgarte_partial_dv.noalias() 
        += coeff * baumgarte_weight_on_velocity_ 
            * frame_a_partial_da_.template topRows<3>();
  }
  if (baumgarte_weight_on_position_ != 0.) {
    baumgarte_partial_dq.noalias()
        += coeff * baumgarte_weight_on_position_ 
            * data.oMf[contact_frame_id_].rotation()
            * J_frame_.template topRows<3>();
  }
}



void PointContact::activate() {
  is_active_ = true;
}


void PointContact::deactivate() {
  is_active_ = false;
}


bool PointContact::isActive() const {
  return is_active_;
}


int PointContact::contact_frame_id() const {
  return contact_frame_id_;
}


int PointContact::parent_joint_id() const {
  return parent_joint_id_;
}


double PointContact::baumgarte_weight_on_velocity() const {
  return baumgarte_weight_on_velocity_;
}


double PointContact::baumgarte_weight_on_position() const {
  return baumgarte_weight_on_position_;
}


Eigen::Vector3d PointContact::contact_point() const {
  return contact_point_;
}

} // namespace idocp