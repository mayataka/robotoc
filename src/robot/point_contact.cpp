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
    fXj_(jXf_.inverse().toActionMatrix().transpose()),
    J_frame_(Eigen::MatrixXd::Zero(6, model.nv)),
    joint_v_partial_dq_(Eigen::MatrixXd::Zero(6, model.nv)),
    joint_a_partial_dq_(Eigen::MatrixXd::Zero(6, model.nv)),
    joint_a_partial_dv_(Eigen::MatrixXd::Zero(6, model.nv)),
    joint_a_partial_da_(Eigen::MatrixXd::Zero(6, model.nv)) {
  assert(contact_frame_id_ >= 0);
  assert(baumgarte_weight_on_velocity_ >= 0);
  assert(baumgarte_weight_on_position_ >= 0);
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
    fXj_(),
    J_frame_(),
    joint_v_partial_dq_(),
    joint_a_partial_dq_(),
    joint_a_partial_dv_(),
    joint_a_partial_da_() {
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
                                      Eigen::MatrixXd& Jacobian, 
                                      const bool transpose) {
  pinocchio::getFrameJacobian(model, data, contact_frame_id_, pinocchio::LOCAL, 
                              J_frame_);
  if (transpose) {
    assert(Jacobian.rows() == dimv_);
    assert(Jacobian.cols() == 3);
    Jacobian = J_frame_.topRows<3>().transpose(); 
  }
  else {
    assert(Jacobian.rows() == 3);
    assert(Jacobian.cols() == dimv_);
    Jacobian = J_frame_.topRows<3>(); 
  }
}


void PointContact::getContactJacobian(const pinocchio::Model& model, 
                                      pinocchio::Data& data, 
                                      const int block_begin_index,
                                      Eigen::MatrixXd& Jacobian,
                                      const bool transpose) {
  assert(block_begin_index >= 0);
  pinocchio::getFrameJacobian(model, data, contact_frame_id_, pinocchio::LOCAL, 
                              J_frame_);
  if (transpose) {
    assert(Jacobian.rows() == dimv_);
    assert(Jacobian.cols() >= 3);
    Jacobian.block(0, block_begin_index, dimv_, 3) 
        = J_frame_.topRows<3>().transpose(); 
  }
  else {
    assert(Jacobian.rows() >= 3);
    assert(Jacobian.cols() == dimv_);
    Jacobian.block(block_begin_index, 0, 3, dimv_) = J_frame_.topRows<3>(); 
  }
}


void PointContact::getContactJacobian(const pinocchio::Model& model, 
                                      pinocchio::Data& data, 
                                      const int block_begin_index, 
                                      const double coeff,
                                      Eigen::MatrixXd& Jacobian,
                                      const bool transpose) {
  assert(block_begin_index >= 0);
  pinocchio::getFrameJacobian(model, data, contact_frame_id_, pinocchio::LOCAL, 
                              J_frame_);
  if (transpose) {
    assert(Jacobian.rows() == dimv_);
    assert(Jacobian.cols() >= 3);
    Jacobian.block(0, block_begin_index, dimv_, 3) 
        = coeff * J_frame_.topRows<3>().transpose(); 
  }
  else {
    assert(Jacobian.rows() >= 3);
    assert(Jacobian.cols() == dimv_);
    Jacobian.block(block_begin_index, 0, 3, dimv_) 
        = coeff * J_frame_.topRows<3>(); 
  }
}


void PointContact::computeBaumgarteResidual(
    const pinocchio::Model& model, const pinocchio::Data& data, 
    Eigen::Vector3d& baumgarte_residual) const {
  assert(baumgarte_residual.size() == 3);
  baumgarte_residual
      = (data.oMi[parent_joint_id_].act(data.a[parent_joint_id_])).linear()
          + baumgarte_weight_on_velocity_
              * (data.oMi[parent_joint_id_].act(data.v[parent_joint_id_])).linear()
          + baumgarte_weight_on_position_
              * (data.oMf[contact_frame_id_].translation()-contact_point_);
}


void PointContact::computeBaumgarteResidual(
    const pinocchio::Model& model, const pinocchio::Data& data, 
    const int result_begin, Eigen::VectorXd& baumgarte_residual) const {
  assert(result_begin >= 0);
  assert(baumgarte_residual.size() >= 3);
  baumgarte_residual.segment<3>(result_begin)
      = (data.oMi[parent_joint_id_].act(data.a[parent_joint_id_])).linear()
          + baumgarte_weight_on_velocity_
              * (data.oMi[parent_joint_id_].act(data.v[parent_joint_id_])).linear()
          + baumgarte_weight_on_position_
              * (data.oMf[contact_frame_id_].translation()-contact_point_);
}


void PointContact::computeBaumgarteResidual(
    const pinocchio::Model& model, const pinocchio::Data& data, 
    const int result_begin, const double coeff,
    Eigen::VectorXd& baumgarte_residual) const {
  assert(result_begin >= 0);
  assert(baumgarte_residual.size() >= 3);
  baumgarte_residual.segment<3>(result_begin)
      = coeff * (data.oMi[parent_joint_id_].act(data.a[parent_joint_id_])).linear()
          + coeff * baumgarte_weight_on_velocity_
              * (data.oMi[parent_joint_id_].act(data.v[parent_joint_id_])).linear()
          + coeff * baumgarte_weight_on_position_
              * (data.oMf[contact_frame_id_].translation()-contact_point_);
}


void PointContact::computeBaumgarteDerivatives(
    const pinocchio::Model& model, pinocchio::Data& data, 
    Eigen::MatrixXd& baumgarte_partial_dq, 
    Eigen::MatrixXd& baumgarte_partial_dv, 
    Eigen::MatrixXd& baumgarte_partial_da) {
  assert(baumgarte_partial_dq.cols() == dimv_);
  assert(baumgarte_partial_dv.cols() == dimv_);
  assert(baumgarte_partial_da.cols() == dimv_);
  assert(baumgarte_partial_dq.rows() == 3);
  assert(baumgarte_partial_dv.rows() == 3);
  assert(baumgarte_partial_da.rows() == 3);
 pinocchio::getJointAccelerationDerivatives(model, data, parent_joint_id_, 
                                            pinocchio::WORLD,
                                            joint_v_partial_dq_, 
                                            joint_a_partial_dq_, 
                                            joint_a_partial_dv_, 
                                            joint_a_partial_da_);
  pinocchio::getFrameJacobian(model, data, contact_frame_id_, 
                              pinocchio::LOCAL_WORLD_ALIGNED, J_frame_);
  baumgarte_partial_dq
      = joint_a_partial_dq_.template topRows<3>()
          + baumgarte_weight_on_velocity_ 
              * joint_v_partial_dq_.template topRows<3>()
          + baumgarte_weight_on_position_ 
              * J_frame_.template topRows<3>();
  baumgarte_partial_dv 
      = joint_a_partial_dv_.template topRows<3>()
          + baumgarte_weight_on_velocity_ 
              * joint_a_partial_da_.template topRows<3>();
  baumgarte_partial_da 
      = joint_a_partial_da_.template topRows<3>();
}


void PointContact::computeBaumgarteDerivatives(
    const pinocchio::Model& model, pinocchio::Data& data, 
    const int block_rows_begin, Eigen::MatrixXd& baumgarte_partial_dq, 
    Eigen::MatrixXd& baumgarte_partial_dv, 
    Eigen::MatrixXd& baumgarte_partial_da) {
  assert(block_rows_begin >= 0);
  assert(baumgarte_partial_dq.cols() == dimv_);
  assert(baumgarte_partial_dv.cols() == dimv_);
  assert(baumgarte_partial_da.cols() == dimv_);
  assert(baumgarte_partial_dq.rows() >= 3);
  assert(baumgarte_partial_dv.rows() >= 3);
  assert(baumgarte_partial_da.rows() >= 3);
  pinocchio::getJointAccelerationDerivatives(model, data, parent_joint_id_, 
                                             pinocchio::WORLD,
                                             joint_v_partial_dq_, 
                                             joint_a_partial_dq_, 
                                             joint_a_partial_dv_, 
                                             joint_a_partial_da_);
  pinocchio::getFrameJacobian(model, data, contact_frame_id_, 
                              pinocchio::LOCAL_WORLD_ALIGNED, J_frame_);
  baumgarte_partial_dq.block(block_rows_begin, 0, 3, dimv_) 
      = joint_a_partial_dq_.template topRows<3>()
          + baumgarte_weight_on_velocity_ 
              * joint_v_partial_dq_.template topRows<3>()
          + baumgarte_weight_on_position_ 
              * J_frame_.template topRows<3>();
  baumgarte_partial_dv.block(block_rows_begin, 0, 3, dimv_) 
      = joint_a_partial_dv_.template topRows<3>()
          + baumgarte_weight_on_velocity_ 
              * joint_a_partial_da_.template topRows<3>();
  baumgarte_partial_da.block(block_rows_begin, 0, 3, dimv_) 
      = joint_a_partial_da_.template topRows<3>();
}


void PointContact::computeBaumgarteDerivatives(
    const pinocchio::Model& model, pinocchio::Data& data, 
    const int block_rows_begin, const double coeff, 
    Eigen::MatrixXd& baumgarte_partial_dq, 
    Eigen::MatrixXd& baumgarte_partial_dv, 
    Eigen::MatrixXd& baumgarte_partial_da) {
  assert(block_rows_begin >= 0);
  assert(baumgarte_partial_dq.cols() == dimv_);
  assert(baumgarte_partial_dv.cols() == dimv_);
  assert(baumgarte_partial_da.cols() == dimv_);
  assert(baumgarte_partial_dq.rows() >= 3);
  assert(baumgarte_partial_dv.rows() >= 3);
  assert(baumgarte_partial_da.rows() >= 3);
  pinocchio::getJointAccelerationDerivatives(model, data, parent_joint_id_, 
                                             pinocchio::WORLD,
                                             joint_v_partial_dq_, 
                                             joint_a_partial_dq_, 
                                             joint_a_partial_dv_, 
                                             joint_a_partial_da_);
  pinocchio::getFrameJacobian(model, data, contact_frame_id_, 
                              pinocchio::LOCAL_WORLD_ALIGNED, J_frame_);
  baumgarte_partial_dq.block(block_rows_begin, 0, 3, dimv_) 
      = coeff * joint_a_partial_dq_.template topRows<3>()
          + coeff * baumgarte_weight_on_velocity_ 
              * joint_v_partial_dq_.template topRows<3>()
          + coeff * baumgarte_weight_on_position_ 
              * J_frame_.template topRows<3>();
  baumgarte_partial_dv.block(block_rows_begin, 0, 3, dimv_) 
      = coeff * joint_a_partial_dv_.template topRows<3>()
          + coeff * baumgarte_weight_on_velocity_ 
              * joint_a_partial_da_.template topRows<3>();
  baumgarte_partial_da.block(block_rows_begin, 0, 3, dimv_) 
      = coeff * joint_a_partial_da_.template topRows<3>();
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


int PointContact::dimv() const {
  return dimv_;
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


pinocchio::SE3 PointContact::jXf() const {
  return jXf_;
}

} // namespace idocp