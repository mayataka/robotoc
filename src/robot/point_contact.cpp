#include "robot/point_contact.hpp"

#include <assert.h>


namespace invdynocp {

PointContact::PointContact(const pinocchio::Model& model, 
                           const unsigned int contact_frame_id, 
                           const double baumgarte_alpha, 
                           const double baumgarte_beta) 
  : contact_frame_id_(contact_frame_id),
    parent_joint_id_(model.frames[contact_frame_id_].parent), 
    dimv_(model.nv),
    baumgarte_alpha_(baumgarte_alpha),
    baumgarte_beta_(baumgarte_beta),
    contact_point_(Eigen::Vector3d::Zero()),
    jXf_(model.frames[contact_frame_id_].placement),
    fXj_(jXf_.inverse().toActionMatrix().transpose()),
    J_frame_(Eigen::MatrixXd::Zero(6, model.nv)),
    joint_v_partial_dq_(Eigen::MatrixXd::Zero(6, model.nv)),
    joint_a_partial_dq_(Eigen::MatrixXd::Zero(6, model.nv)),
    joint_a_partial_dv_(Eigen::MatrixXd::Zero(6, model.nv)),
    joint_a_partial_da_(Eigen::MatrixXd::Zero(6, model.nv)) {
}


PointContact::~PointContact() {
}


PointContact::PointContact(PointContact&& other) noexcept
  : contact_frame_id_(other.contact_frame_id()),
    parent_joint_id_(other.parent_joint_id()), 
    dimv_(other.dimv()),
    baumgarte_alpha_(other.baumgarte_alpha()),
    baumgarte_beta_(other.baumgarte_beta()),
    contact_point_(other.contact_point()),
    jXf_(other.jXf()),
    fXj_(other.jXf().inverse().toActionMatrix().transpose()),
    J_frame_(),
    joint_v_partial_dq_(),
    joint_a_partial_dq_(),
    joint_a_partial_dv_(),
    joint_a_partial_da_() {
  J_frame_.resize(6, dimv_); 
  joint_v_partial_dq_.resize(6, dimv_); 
  joint_a_partial_dq_.resize(6, dimv_); 
  joint_a_partial_dv_.resize(6, dimv_); 
  joint_a_partial_da_.resize(6, dimv_); 
}


PointContact& PointContact::operator=(PointContact&& other) noexcept {
  contact_frame_id_ = other.contact_frame_id();
  parent_joint_id_ = other.parent_joint_id();
  dimv_ = other.dimv();
  baumgarte_alpha_ = other.baumgarte_alpha();
  baumgarte_beta_ = other.baumgarte_beta();
  contact_point_ = other.contact_point();
  jXf_ = other.jXf();
  fXj_ = other.jXf().inverse().toActionMatrix().transpose();
  J_frame_.resize(6, dimv_); 
  joint_v_partial_dq_.resize(6, dimv_); 
  joint_a_partial_dq_.resize(6, dimv_); 
  joint_a_partial_dv_.resize(6, dimv_); 
  joint_a_partial_da_.resize(6, dimv_); 
}


void PointContact::resetBaugrarteParameters(const double alpha, 
                                            const double beta) {
  baumgarte_alpha_ = alpha;
  baumgarte_beta_ = beta;
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
                                      Eigen::MatrixXd& J_contact) const {
  pinocchio::getFrameJacobian(model, data, contact_frame_id_, pinocchio::LOCAL, 
                              J_frame_);
  J_contact = J_frame_.topRows<3>(); 
}


void PointContact::getContactJacobian(const pinocchio::Model& model, 
                                      pinocchio::Data& data, 
                                      const unsigned int result_mat_row_begin, 
                                      const unsigned int result_mat_col_begin, 
                                      Eigen::MatrixXd& J_contacts) const {
  pinocchio::getFrameJacobian(model, data, contact_frame_id_, pinocchio::LOCAL, 
                              J_frame_);
  J_contacts.block(result_mat_row_begin, result_mat_col_begin, 3, dimv_) 
      = J_frame_.topRows<3>(); 
}


void PointContact::computeBaumgarteResidual(
    const pinocchio::Model& model, const pinocchio::Data& data, 
    Eigen::Vector3d& baumgarte_residual) const {
  baumgarte_residual
      = (data.oMi[parent_joint_id_].act(data.a[parent_joint_id_])).linear()
          + baumgarte_alpha_ 
              * (data.oMi[parent_joint_id_].act(data.v[parent_joint_id_])).linear()
          + baumgarte_beta_  
              * (data.oMf[contact_frame_id_].translation()-contact_point_);
}


void PointContact::computeBaumgarteResidual(
    const pinocchio::Model& model, const pinocchio::Data& data, 
    const unsigned int result_vec_start_index, 
    Eigen::VectorXd& baumgarte_residual) const {
  baumgarte_residual.segment<3>(result_vec_start_index)
      = (data.oMi[parent_joint_id_].act(data.a[parent_joint_id_])).linear()
          + baumgarte_alpha_ 
              * (data.oMi[parent_joint_id_].act(data.v[parent_joint_id_])).linear()
          + baumgarte_beta_  
              * (data.oMf[contact_frame_id_].translation()-contact_point_);
}


void PointContact::computeBaumgarteDerivatives(
    const pinocchio::Model& model, pinocchio::Data& data, 
    Eigen::MatrixXd& baumgarte_partial_dq, 
    Eigen::MatrixXd& baumgarte_partial_dv, 
    Eigen::MatrixXd& baumgarte_partial_da) {
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
          + baumgarte_alpha_ * joint_v_partial_dq_.template topRows<3>()
          + baumgarte_beta_ * J_frame_.template topRows<3>();
  baumgarte_partial_dv 
      = joint_a_partial_dv_.template topRows<3>()
          + baumgarte_alpha_ * joint_a_partial_da_.template topRows<3>();
  baumgarte_partial_da 
      = joint_a_partial_da_.template topRows<3>();
}


void PointContact::computeBaumgarteDerivatives(
    const pinocchio::Model& model, pinocchio::Data& data, 
    const unsigned int result_mat_row_begin, 
    const unsigned int result_mat_col_begin,
    Eigen::MatrixXd& baumgarte_partial_dq, 
    Eigen::MatrixXd& baumgarte_partial_dv, 
    Eigen::MatrixXd& baumgarte_partial_da) {
  pinocchio::getJointAccelerationDerivatives(model, data, parent_joint_id_, 
                                             pinocchio::WORLD,
                                             joint_v_partial_dq_, 
                                             joint_a_partial_dq_, 
                                             joint_a_partial_dv_, 
                                             joint_a_partial_da_);
  pinocchio::getFrameJacobian(model, data, contact_frame_id_, 
                              pinocchio::LOCAL_WORLD_ALIGNED, J_frame_);
  baumgarte_partial_dq.block(result_mat_row_begin, result_mat_col_begin, 3, 
                             dimv_) 
      = joint_a_partial_dq_.template topRows<3>()
          + baumgarte_alpha_ * joint_v_partial_dq_.template topRows<3>()
          + baumgarte_beta_ * J_frame_.template topRows<3>();
  baumgarte_partial_dv.block(result_mat_row_begin, result_mat_col_begin, 3, 
                             dimv_) 
      = joint_a_partial_dv_.template topRows<3>()
          + baumgarte_alpha_ * joint_a_partial_da_.template topRows<3>();
  baumgarte_partial_da.block(result_mat_row_begin, result_mat_col_begin, 3, 
                             dimv_) 
      = joint_a_partial_da_.template topRows<3>();
}


unsigned int PointContact::contact_frame_id() const {
  return contact_frame_id_;
}


unsigned int PointContact::parent_joint_id() const {
  return parent_joint_id_;
}


unsigned int PointContact::dimv() const {
  return dimv_;
}


double PointContact::baumgarte_alpha() const {
  return baumgarte_alpha_;
}


double PointContact::baumgarte_beta() const {
  return baumgarte_beta_;
}


Eigen::Vector3d PointContact::contact_point() const {
  return contact_point_;
}


pinocchio::SE3 PointContact::jXf() const {
  return jXf_;
}

} // namespace invdynocp 