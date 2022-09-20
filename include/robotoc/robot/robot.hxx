#ifndef ROBOTOC_ROBOT_HXX_
#define ROBOTOC_ROBOT_HXX_

#include "robotoc/robot/robot.hpp"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/kinematics-derivatives.hpp"
#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/frames-derivatives.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/rnea-derivatives.hpp"
#include "pinocchio/algorithm/aba.hpp"
#include "pinocchio/algorithm/cholesky.hpp"
#include "pinocchio/algorithm/contact-dynamics.hpp"

#include <cassert>
#include <stdexcept>

namespace robotoc {

template <typename TangentVectorType, typename ConfigVectorType>
inline void Robot::integrateConfiguration(
    const Eigen::MatrixBase<TangentVectorType>& v, 
    const double integration_length, 
    const Eigen::MatrixBase<ConfigVectorType>& q) const {
  assert(v.size() == dimv_);
  assert(q.size() == dimq_);
  if (info_.base_joint_type == BaseJointType::FloatingBase) {
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
  pinocchio::integrate(
      model_, q, integration_length*v, 
      const_cast<Eigen::MatrixBase<ConfigVectorType2>&>(q_integrated));
}


template <typename ConfigVectorType, typename TangentVectorType, 
          typename MatrixType1, typename MatrixType2>
inline void Robot::dIntegrateTransport_dq(
    const Eigen::MatrixBase<ConfigVectorType>& q,
    const Eigen::MatrixBase<TangentVectorType>& v,
    const Eigen::MatrixBase<MatrixType1>& Jin,
    const Eigen::MatrixBase<MatrixType2>& Jout) const {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(Jin.cols() == dimv_);
  assert(Jout.rows() == Jin.rows());
  assert(Jout.cols() == Jin.cols());
  pinocchio::dIntegrateTransport(
      model_, q, v, Jin.transpose(), 
      const_cast<Eigen::MatrixBase<MatrixType2>&>(Jout).transpose(),
      pinocchio::ARG0);
}


template <typename ConfigVectorType, typename TangentVectorType, 
          typename MatrixType1, typename MatrixType2>
inline void Robot::dIntegrateTransport_dv(
    const Eigen::MatrixBase<ConfigVectorType>& q,
    const Eigen::MatrixBase<TangentVectorType>& v,
    const Eigen::MatrixBase<MatrixType1>& Jin,
    const Eigen::MatrixBase<MatrixType2>& Jout) const {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(Jin.cols() == dimv_);
  assert(Jout.rows() == Jin.rows());
  assert(Jout.cols() == Jin.cols());
  pinocchio::dIntegrateTransport(
      model_, q, v, Jin.transpose(), 
      const_cast<Eigen::MatrixBase<MatrixType2>&>(Jout).transpose(),
      pinocchio::ARG1);
}


template <typename ConfigVectorType1, typename ConfigVectorType2, 
          typename TangentVectorType>
inline void Robot::subtractConfiguration(
    const Eigen::MatrixBase<ConfigVectorType1>& qf, 
    const Eigen::MatrixBase<ConfigVectorType2>& q0,
    const Eigen::MatrixBase<TangentVectorType>& qdiff) const {
  assert(qf.size() == dimq_);
  assert(q0.size() == dimq_);
  assert(qdiff.size() == dimv_);
  pinocchio::difference(
      model_, q0, qf, 
      const_cast<Eigen::MatrixBase<TangentVectorType>&>(qdiff));
}


template <typename ConfigVectorType1, typename ConfigVectorType2, 
          typename MatrixType>
inline void Robot::dSubtractConfiguration_dqf(
    const Eigen::MatrixBase<ConfigVectorType1>& qf,
    const Eigen::MatrixBase<ConfigVectorType2>& q0,
    const Eigen::MatrixBase<MatrixType>& dqdiff_dqf) const {
  assert(qf.size() == dimq_);
  assert(q0.size() == dimq_);
  assert(dqdiff_dqf.rows() == dimv_);
  assert(dqdiff_dqf.cols() == dimv_);
  pinocchio::dDifference(model_, q0, qf, 
                         const_cast<Eigen::MatrixBase<MatrixType>&>(dqdiff_dqf),
                         pinocchio::ARG1);
}


template <typename ConfigVectorType1, typename ConfigVectorType2, 
          typename MatrixType>
inline void Robot::dSubtractConfiguration_dq0(
    const Eigen::MatrixBase<ConfigVectorType1>& qf,
    const Eigen::MatrixBase<ConfigVectorType2>& q0,
    const Eigen::MatrixBase<MatrixType>& dqdiff_dq0) const {
  assert(qf.size() == dimq_);
  assert(q0.size() == dimq_);
  assert(dqdiff_dq0.rows() == dimv_);
  assert(dqdiff_dq0.cols() == dimv_);
  pinocchio::dDifference(model_, q0, qf, 
                         const_cast<Eigen::MatrixBase<MatrixType>&>(dqdiff_dq0),
                         pinocchio::ARG0);
}


template <typename ConfigVectorType1, typename ConfigVectorType2, 
          typename ConfigVectorType3>
inline void Robot::interpolateConfiguration(
    const Eigen::MatrixBase<ConfigVectorType1>& q1, 
    const Eigen::MatrixBase<ConfigVectorType2>& q2, const double t,
    const Eigen::MatrixBase<ConfigVectorType3>& qout) const {
  assert(q1.size() == dimq_);
  assert(q2.size() == dimq_);
  assert(qout.size() == dimq_);
  assert(t >= 0.0);
  assert(t <= 1.0);
  pinocchio::interpolate(model_, q1, q2, t, 
                         const_cast<Eigen::MatrixBase<ConfigVectorType3>&>(qout));
}


template <typename ConfigVectorType, typename MatrixType>
inline void Robot::integrateCoeffWiseJacobian(
    const Eigen::MatrixBase<ConfigVectorType>& q, 
    const Eigen::MatrixBase<MatrixType>& J) const {
  assert(q.size() == dimq_);
  assert(J.rows() == dimq_);
  assert(J.cols() == dimv_);
  pinocchio::integrateCoeffWiseJacobian(model_, q, 
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
  pinocchio::jacobianCenterOfMass(model_, data_, false);
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
  pinocchio::jacobianCenterOfMass(model_, data_, false);
}


template <typename ConfigVectorType>
inline void Robot::updateKinematics(
    const Eigen::MatrixBase<ConfigVectorType>& q) {
  assert(q.size() == dimq_);
  pinocchio::framesForwardKinematics(model_, data_, q);
  pinocchio::computeJointJacobians(model_, data_, q);
  pinocchio::jacobianCenterOfMass(model_, data_, false);
}


template <typename ConfigVectorType, typename TangentVectorType1, 
          typename TangentVectorType2>
inline void Robot::updateFrameKinematics(
    const Eigen::MatrixBase<ConfigVectorType>& q, 
    const Eigen::MatrixBase<TangentVectorType1>& v, 
    const Eigen::MatrixBase<TangentVectorType2>& a) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  pinocchio::forwardKinematics(model_, data_, q, v, a);
  pinocchio::updateFramePlacements(model_, data_);
  pinocchio::centerOfMass(model_, data_, q, v, a, false);
}


template <typename ConfigVectorType, typename TangentVectorType>
inline void Robot::updateFrameKinematics(
    const Eigen::MatrixBase<ConfigVectorType>& q, 
    const Eigen::MatrixBase<TangentVectorType>& v) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  pinocchio::forwardKinematics(model_, data_, q, v);
  pinocchio::updateFramePlacements(model_, data_);
  pinocchio::centerOfMass(model_, data_, q, v, false);
}


template <typename ConfigVectorType>
inline void Robot::updateFrameKinematics(
    const Eigen::MatrixBase<ConfigVectorType>& q) {
  assert(q.size() == dimq_);
  pinocchio::framesForwardKinematics(model_, data_, q);
  pinocchio::centerOfMass(model_, data_, q, false);
}


inline const Eigen::Vector3d& Robot::framePosition(const int frame_id) const {
  return data_.oMf[frame_id].translation();
}


inline const Eigen::Vector3d& Robot::framePosition(
    const std::string& frame_name) const {
  return framePosition(frameId(frame_name));
}


inline const Eigen::Matrix3d& Robot::frameRotation(const int frame_id) const {
  return data_.oMf[frame_id].rotation();
}


inline const Eigen::Matrix3d& Robot::frameRotation(
    const std::string& frame_name) const {
  return frameRotation(frameId(frame_name));
}


inline const SE3& Robot::framePlacement(const int frame_id) const {
  return data_.oMf[frame_id];
}


inline const SE3& Robot::framePlacement(const std::string& frame_name) const {
  return framePlacement(frameId(frame_name));
}


inline const Eigen::Vector3d& Robot::CoM() const {
  return data_.com[0];
}


inline Eigen::Vector3d Robot::frameLinearVelocity(
    const int frame_id, const pinocchio::ReferenceFrame reference_frame) const {
  return pinocchio::getFrameVelocity(model_, data_, frame_id, reference_frame).linear();
}


inline Eigen::Vector3d Robot::frameLinearVelocity(
    const std::string& frame_name, 
    const pinocchio::ReferenceFrame reference_frame) const {
  return frameLinearVelocity(frameId(frame_name), reference_frame);
}


inline Eigen::Vector3d Robot::frameAngularVelocity(
    const int frame_id, const pinocchio::ReferenceFrame reference_frame) const {
  return pinocchio::getFrameVelocity(model_, data_, frame_id, reference_frame).angular();
}


inline Eigen::Vector3d Robot::frameAngularVelocity(
    const std::string& frame_name, 
    const pinocchio::ReferenceFrame reference_frame) const {
  return frameAngularVelocity(frameId(frame_name), reference_frame);
}


inline Robot::Vector6d Robot::frameSpatialVelocity(
    const int frame_id, const pinocchio::ReferenceFrame reference_frame) const {
  return pinocchio::getFrameVelocity(model_, data_, frame_id, reference_frame).toVector();
}


inline Robot::Vector6d Robot::frameSpatialVelocity(
    const std::string& frame_name, 
    const pinocchio::ReferenceFrame reference_frame) const {
  return frameSpatialVelocity(frameId(frame_name), reference_frame);
}


inline const Eigen::Vector3d& Robot::CoMVelocity() const {
  return data_.vcom[0];
}


template <typename MatrixType>
inline void Robot::getFrameJacobian(const int frame_id, 
                                    const Eigen::MatrixBase<MatrixType>& J) {
  assert(J.rows() == 6);
  assert(J.cols() == dimv_);
  pinocchio::getFrameJacobian(model_, data_, frame_id, pinocchio::LOCAL, 
                              const_cast<Eigen::MatrixBase<MatrixType>&>(J));
}


template <typename MatrixType>
inline void Robot::getCoMJacobian(const Eigen::MatrixBase<MatrixType>& J) const {
  assert(J.rows() == 3);
  assert(J.cols() == dimv_);
  const_cast<Eigen::MatrixBase<MatrixType>&>(J) = data_.Jcom;
}


template <typename Vector3dType>
inline void Robot::transformFromLocalToWorld(
    const int frame_id, const Eigen::Vector3d& vec_local, 
    const Eigen::MatrixBase<Vector3dType>& vec_world) const {
  assert(vec_world.size() == 3);
  const_cast<Eigen::MatrixBase<Vector3dType>&>(vec_world).noalias()
      = frameRotation(frame_id) * vec_local;
}


template <typename Vector3dType, typename MatrixType>
inline void Robot::getJacobianTransformFromLocalToWorld(
    const int frame_id, const Eigen::MatrixBase<Vector3dType>& vec_world,
    const Eigen::MatrixBase<MatrixType>& J) {
  assert(vec_world.size() == 3);
  assert(J.rows() == 6);
  assert(J.cols() == dimv_);
  const_cast<Eigen::MatrixBase<MatrixType>&>(J).setZero();
  getFrameJacobian(frame_id, const_cast<Eigen::MatrixBase<MatrixType>&>(J));
  for (int i=0; i<dimv_; ++i) {
    const_cast<Eigen::MatrixBase<MatrixType>&>(J).template topRows<3>().col(i).noalias()
        = J.template bottomRows<3>().col(i).cross(vec_world.template head<3>());
  }
}


template <typename VectorType>
inline void Robot::computeBaumgarteResidual(
    const ContactStatus& contact_status, 
    const Eigen::MatrixBase<VectorType>& baumgarte_residual) {
  assert(baumgarte_residual.size() == contact_status.dimf());
  const int num_point_contacts = point_contacts_.size();
  const int num_surface_contacts = surface_contacts_.size();
  int dimf = 0;
  for (int i=0; i<num_point_contacts; ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_[i].computeBaumgarteResidual(
          model_, data_, contact_status.contactPosition(i),
          (const_cast<Eigen::MatrixBase<VectorType>&>(baumgarte_residual))
              .template segment<3>(dimf));
      dimf += 3;
    }
  }
  for (int i=0; i<num_surface_contacts; ++i) {
    if (contact_status.isContactActive(i+num_point_contacts)) {
      surface_contacts_[i].computeBaumgarteResidual(
          model_, data_, contact_status.contactPlacement(i+num_point_contacts),
          (const_cast<Eigen::MatrixBase<VectorType>&>(baumgarte_residual))
              .template segment<6>(dimf));
      dimf += 6;
    }
  }
  assert(dimf == contact_status.dimf());
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void Robot::computeBaumgarteDerivatives(
    const ContactStatus& contact_status,
    const Eigen::MatrixBase<MatrixType1>& baumgarte_partial_dq, 
    const Eigen::MatrixBase<MatrixType2>& baumgarte_partial_dv, 
    const Eigen::MatrixBase<MatrixType3>& baumgarte_partial_da) {
  assert(baumgarte_partial_dq.rows() == contact_status.dimf());
  assert(baumgarte_partial_dq.cols() == dimv_);
  assert(baumgarte_partial_dv.rows() == contact_status.dimf());
  assert(baumgarte_partial_dv.cols() == dimv_);
  assert(baumgarte_partial_da.rows() == contact_status.dimf());
  assert(baumgarte_partial_da.cols() == dimv_);
  const int num_point_contacts = point_contacts_.size();
  const int num_surface_contacts = surface_contacts_.size();
  int dimf = 0;
  for (int i=0; i<num_point_contacts; ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_[i].computeBaumgarteDerivatives(
          model_, data_, 
          (const_cast<Eigen::MatrixBase<MatrixType1>&>(baumgarte_partial_dq))
              .block(dimf, 0, 3, dimv_),
          (const_cast<Eigen::MatrixBase<MatrixType2>&>(baumgarte_partial_dv))
              .block(dimf, 0, 3, dimv_),
          (const_cast<Eigen::MatrixBase<MatrixType3>&>(baumgarte_partial_da))
              .block(dimf, 0, 3, dimv_));
      dimf += 3;
    }
  }
  for (int i=0; i<num_surface_contacts; ++i) {
    if (contact_status.isContactActive(i+num_point_contacts)) {
      surface_contacts_[i].computeBaumgarteDerivatives(
          model_, data_, 
          (const_cast<Eigen::MatrixBase<MatrixType1>&>(baumgarte_partial_dq))
              .block(dimf, 0, 6, dimv_),
          (const_cast<Eigen::MatrixBase<MatrixType2>&>(baumgarte_partial_dv))
              .block(dimf, 0, 6, dimv_),
          (const_cast<Eigen::MatrixBase<MatrixType3>&>(baumgarte_partial_da))
              .block(dimf, 0, 6, dimv_));
      dimf += 6;
    }
  }
  assert(dimf == contact_status.dimf());
}


template <typename VectorType>
inline void Robot::computeImpulseVelocityResidual(
    const ImpulseStatus& impulse_status, 
    const Eigen::MatrixBase<VectorType>& velocity_residual) const {
  assert(velocity_residual.size() == impulse_status.dimi());
  const int num_point_contacts = point_contacts_.size();
  const int num_surface_contacts = surface_contacts_.size();
  int dimf = 0;
  for (int i=0; i<num_point_contacts; ++i) {
    if (impulse_status.isImpulseActive(i)) {
      point_contacts_[i].computeContactVelocityResidual(
          model_, data_, 
          (const_cast<Eigen::MatrixBase<VectorType>&>(velocity_residual))
              .template segment<3>(dimf));
      dimf += 3;
    }
  }
  for (int i=0; i<num_surface_contacts; ++i) {
    if (impulse_status.isImpulseActive(i+num_point_contacts)) {
      surface_contacts_[i].computeContactVelocityResidual(
          model_, data_, 
          (const_cast<Eigen::MatrixBase<VectorType>&>(velocity_residual))
              .template segment<6>(dimf));
      dimf += 6;
    }
  }
  assert(dimf == impulse_status.dimi());
}


template <typename MatrixType1, typename MatrixType2>
inline void Robot::computeImpulseVelocityDerivatives(
    const ImpulseStatus& impulse_status, 
    const Eigen::MatrixBase<MatrixType1>& velocity_partial_dq, 
    const Eigen::MatrixBase<MatrixType2>& velocity_partial_dv) {
  assert(velocity_partial_dq.rows() == impulse_status.dimi());
  assert(velocity_partial_dq.cols() == dimv_);
  assert(velocity_partial_dv.rows() == impulse_status.dimi());
  assert(velocity_partial_dv.cols() == dimv_);
  const int num_point_contacts = point_contacts_.size();
  const int num_surface_contacts = surface_contacts_.size();
  int dimf = 0;
  for (int i=0; i<num_point_contacts; ++i) {
    if (impulse_status.isImpulseActive(i)) {
      point_contacts_[i].computeContactVelocityDerivatives(
          model_, data_, 
          (const_cast<Eigen::MatrixBase<MatrixType1>&>(velocity_partial_dq))
              .block(dimf, 0, 3, dimv_),
          (const_cast<Eigen::MatrixBase<MatrixType2>&>(velocity_partial_dv))
              .block(dimf, 0, 3, dimv_));
      dimf += 3;
    }
  }
  for (int i=0; i<num_surface_contacts; ++i) {
    if (impulse_status.isImpulseActive(i+num_point_contacts)) {
      surface_contacts_[i].computeContactVelocityDerivatives(
          model_, data_, 
          (const_cast<Eigen::MatrixBase<MatrixType1>&>(velocity_partial_dq))
              .block(dimf, 0, 6, dimv_),
          (const_cast<Eigen::MatrixBase<MatrixType2>&>(velocity_partial_dv))
              .block(dimf, 0, 6, dimv_));
      dimf += 6;
    }
  }
  assert(dimf == impulse_status.dimi());
}


template <typename VectorType>
inline void Robot::computeContactPositionResidual(
    const ImpulseStatus& impulse_status, 
    const Eigen::MatrixBase<VectorType>& position_residual) {
  assert(position_residual.size() == impulse_status.dimi());
  const int num_point_contacts = point_contacts_.size();
  const int num_surface_contacts = surface_contacts_.size();
  int dimf = 0;
  for (int i=0; i<num_point_contacts; ++i) {
    if (impulse_status.isImpulseActive(i)) {
      point_contacts_[i].computeContactPositionResidual(
          model_, data_, impulse_status.contactPosition(i),
          (const_cast<Eigen::MatrixBase<VectorType>&>(position_residual))
              .template segment<3>(dimf));
      dimf += 3;
    }
  }
  for (int i=0; i<num_surface_contacts; ++i) {
    if (impulse_status.isImpulseActive(i+num_point_contacts)) {
      surface_contacts_[i].computeContactPositionResidual(
          model_, data_, impulse_status.contactPlacement(i+num_point_contacts),
          (const_cast<Eigen::MatrixBase<VectorType>&>(position_residual))
              .template segment<6>(dimf));
      dimf += 6;
    }
  }
  assert(dimf == impulse_status.dimi());
}


template <typename MatrixType>
inline void Robot::computeContactPositionDerivative(
    const ImpulseStatus& impulse_status, 
    const Eigen::MatrixBase<MatrixType>& position_partial_dq) {
  assert(position_partial_dq.rows() == impulse_status.dimi());
  assert(position_partial_dq.cols() == dimv_);
  const int num_point_contacts = point_contacts_.size();
  const int num_surface_contacts = surface_contacts_.size();
  int dimf = 0;
  for (int i=0; i<num_point_contacts; ++i) {
    if (impulse_status.isImpulseActive(i)) {
      point_contacts_[i].computeContactPositionDerivative(
          model_, data_, 
        (const_cast<Eigen::MatrixBase<MatrixType>&>(position_partial_dq))
            .block(dimf, 0, 3, dimv_));
      dimf += 3;
    }
  }
  for (int i=0; i<num_surface_contacts; ++i) {
    if (impulse_status.isImpulseActive(i+num_point_contacts)) {
      surface_contacts_[i].computeContactPositionDerivative(
          model_, data_, 
        (const_cast<Eigen::MatrixBase<MatrixType>&>(position_partial_dq))
            .block(dimf, 0, 6, dimv_));
      dimf += 6;
    }
  }
  assert(dimf == impulse_status.dimi());
}


inline void Robot::setContactForces(const ContactStatus& contact_status, 
                                    const std::vector<Vector6d>& f) {
  assert(f.size() == max_num_contacts_);
  const int num_point_contacts = point_contacts_.size();
  const int num_surface_contacts = surface_contacts_.size();
  for (int i=0; i<num_point_contacts; ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_[i].computeJointForceFromContactForce(
          f[i].template head<3>(), fjoint_);
    }
    else {
      point_contacts_[i].computeJointForceFromContactForce(
          Eigen::Vector3d::Zero(), fjoint_);
    }
  }
  for (int i=0; i<num_surface_contacts; ++i) {
    if (contact_status.isContactActive(i+num_point_contacts)) {
      surface_contacts_[i].computeJointForceFromContactWrench(
          f[i+num_point_contacts], fjoint_);
    }
    else {
      surface_contacts_[i].computeJointForceFromContactWrench(
          Vector6d::Zero(), fjoint_);
    }
  }
}


inline void Robot::setImpulseForces(const ImpulseStatus& impulse_status, 
                                    const std::vector<Vector6d>& f) {
  assert(f.size() == max_num_contacts_);
  const int num_point_contacts = point_contacts_.size();
  const int num_surface_contacts = surface_contacts_.size();
  for (int i=0; i<num_point_contacts; ++i) {
    if (impulse_status.isImpulseActive(i)) {
      point_contacts_[i].computeJointForceFromContactForce(
          f[i].template head<3>(), fjoint_);
    }
    else {
      point_contacts_[i].computeJointForceFromContactForce(
          Eigen::Vector3d::Zero(), fjoint_);
    }
  }
  for (int i=0; i<num_surface_contacts; ++i) {
    if (impulse_status.isImpulseActive(i+num_point_contacts)) {
      surface_contacts_[i].computeJointForceFromContactWrench(
          f[i+num_point_contacts], fjoint_);
    }
    else {
      surface_contacts_[i].computeJointForceFromContactWrench(
          Vector6d::Zero(), fjoint_);
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
  if (max_num_contacts_) {
    const_cast<Eigen::MatrixBase<TangentVectorType3>&>(tau)
        = pinocchio::rnea(model_, data_, q, v, a, fjoint_);
  }
  else {
    const_cast<Eigen::MatrixBase<TangentVectorType3>&>(tau)
        = pinocchio::rnea(model_, data_, q, v, a);
  }
  if (properties_.has_generalized_momentum_bias) {
    const_cast<Eigen::MatrixBase<TangentVectorType3>&>(tau).noalias()
        -= properties_.generalized_momentum_bias;
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
  if (max_num_contacts_) {
    pinocchio::computeRNEADerivatives(
        model_, data_, q, v, a, fjoint_,
        const_cast<Eigen::MatrixBase<MatrixType1>&>(dRNEA_partial_dq),
        const_cast<Eigen::MatrixBase<MatrixType2>&>(dRNEA_partial_dv),
        const_cast<Eigen::MatrixBase<MatrixType3>&>(dRNEA_partial_da));
  }
  else {
    pinocchio::computeRNEADerivatives(
        model_, data_, q, v, a, 
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


template <typename MatrixType1, typename MatrixType2>
inline void Robot::computeMinv(const Eigen::MatrixBase<MatrixType1>& M, 
                               const Eigen::MatrixBase<MatrixType2>& Minv) {
  assert(M.rows() == dimv_);
  assert(M.cols() == dimv_);
  assert(Minv.rows() == dimv_);
  assert(Minv.cols() == dimv_);
  data_.M = M;
  pinocchio::cholesky::decompose(model_, data_);
  pinocchio::cholesky::computeMinv(
      model_, data_, const_cast<Eigen::MatrixBase<MatrixType2>&>(Minv));
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void Robot::computeMJtJinv(
    const Eigen::MatrixBase<MatrixType1>& M, 
    const Eigen::MatrixBase<MatrixType2>& J, 
    const Eigen::MatrixBase<MatrixType3>& MJtJinv) {
  assert(M.rows() == dimv_);
  assert(M.cols() == dimv_);
  assert(J.rows() <= max_dimf_);
  assert(J.cols() == dimv_);
  assert(MJtJinv.rows() == M.rows()+J.rows());
  assert(MJtJinv.cols() == M.rows()+J.rows());
  const int dimf = J.rows();
  data_.M = M;
  pinocchio::cholesky::decompose(model_, data_);
  data_.sDUiJt.leftCols(dimf) = J.transpose();
  pinocchio::cholesky::Uiv(model_, data_, data_.sDUiJt.leftCols(dimf));
  for (Eigen::DenseIndex k=0; k<dimv_; ++k) {
    data_.sDUiJt.leftCols(dimf).row(k) /= std::sqrt(data_.D[k]);
  }
  data_.JMinvJt.topLeftCorner(dimf, dimf).noalias() 
      = data_.sDUiJt.leftCols(dimf).transpose() * data_.sDUiJt.leftCols(dimf);
  if (info_.contact_inv_damping > 0.) {
    data_.JMinvJt.diagonal().array() += info_.contact_inv_damping;
  }
  data_.llt_JMinvJt.compute(data_.JMinvJt.topLeftCorner(dimf, dimf));
  assert(data_.llt_JMinvJt.info() == Eigen::Success);
  Eigen::Block<MatrixType3> topLeft 
      = const_cast<Eigen::MatrixBase<MatrixType3>&>(MJtJinv).topLeftCorner(dimv_, dimv_);
  Eigen::Block<MatrixType3> topRight 
      = const_cast<Eigen::MatrixBase<MatrixType3>&>(MJtJinv).topRightCorner(dimv_, dimf);
  Eigen::Block<MatrixType3> bottomLeft 
      = const_cast<Eigen::MatrixBase<MatrixType3>&>(MJtJinv).bottomLeftCorner(dimf, dimv_);
  Eigen::Block<MatrixType3> bottomRight 
      = const_cast<Eigen::MatrixBase<MatrixType3>&>(MJtJinv).bottomRightCorner(dimf, dimf);
  bottomRight = - pinocchio::Data::MatrixXs::Identity(dimf, dimf);
  topLeft.setIdentity();
  data_.llt_JMinvJt.solveInPlace(bottomRight);
  pinocchio::cholesky::solve(model_, data_, topLeft);
  bottomLeft.noalias() = J * topLeft;
  topRight.noalias() = bottomLeft.transpose() * (-bottomRight);
  topLeft.noalias() -= topRight*bottomLeft;
  bottomLeft = topRight.transpose();
  assert(!MJtJinv.hasNaN());
}


inline Eigen::VectorXd Robot::generateFeasibleConfiguration() const {
  Eigen::VectorXd q_min(dimq_), q_max(dimq_);
  if (info_.base_joint_type == BaseJointType::FloatingBase) {
    q_min.template head<7>() = - Eigen::VectorXd::Ones(7);
    q_max.template head<7>() = Eigen::VectorXd::Ones(7);
  }
  q_min.tail(dimu_) = lower_joint_position_limit_;
  q_max.tail(dimu_) = upper_joint_position_limit_;
  return pinocchio::randomConfiguration(model_, q_min, q_max);
}


template <typename ConfigVectorType>
inline void Robot::normalizeConfiguration(
    const Eigen::MatrixBase<ConfigVectorType>& q) const {
  assert(q.size() == dimq_);
  if (info_.base_joint_type == BaseJointType::FloatingBase) {
    if (q.template segment<4>(3).squaredNorm() 
          <= std::numeric_limits<double>::epsilon()) {
      (const_cast<Eigen::MatrixBase<ConfigVectorType>&> (q)).coeffRef(3) = 1;
    }
    pinocchio::normalize(model_, 
                         const_cast<Eigen::MatrixBase<ConfigVectorType>&>(q));
  }
}


inline int Robot::frameId(const std::string& frame_name) const {
  if (!model_.existFrame(frame_name)) {
    throw std::invalid_argument(
        "[Robot] invalid argument: frame '" + frame_name + "' does not exit!");
  }
  return model_.getFrameId(frame_name);
}


inline std::string Robot::frameName(const int frame_id) const {
  return  model_.frames[frame_id].name;
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


inline double Robot::totalMass() const {
  return pinocchio::computeTotalMass(model_);
}


inline double Robot::totalWeight() const {
  return (- pinocchio::computeTotalMass(model_) * model_.gravity981.coeff(2));
}


inline int Robot::dimq() const {
  return dimq_;
}


inline int Robot::dimv() const {
  return dimv_;
}


inline int Robot::dimu() const {
  return dimu_;
}


inline int Robot::max_dimf() const {
  return max_dimf_;
}


inline int Robot::dim_passive() const {
  return dim_passive_;
}


inline bool Robot::hasFloatingBase() const {
  return (info_.base_joint_type == BaseJointType::FloatingBase);
}


inline int Robot::maxNumContacts() const {
  return (maxNumPointContacts() + maxNumSurfaceContacts());
}


inline int Robot::maxNumPointContacts() const {
  return point_contacts_.size();
}


inline int Robot::maxNumSurfaceContacts() const {
  return surface_contacts_.size();
}


inline ContactType Robot::contactType(const int contact_index) const {
  assert(contact_index >= 0);
  assert(contact_index < max_num_contacts_);
  if (contact_index < maxNumPointContacts()) {
    return ContactType::PointContact;
  }
  else {
    return ContactType::SurfaceContact;
  }
}


inline std::vector<ContactType> Robot::contactTypes() const {
  std::vector<ContactType> contact_types;
  for (const auto& e : point_contacts_) {
    contact_types.push_back(ContactType::PointContact);
  }
  for (const auto& e : surface_contacts_) {
    contact_types.push_back(ContactType::SurfaceContact);
  }
  return contact_types;
}


inline std::vector<int> Robot::contactFrames() const {
  std::vector<int> contact_frames;
  for (const auto& e : point_contacts_) {
    contact_frames.push_back(e.contactFrameId());
  }
  for (const auto& e : surface_contacts_) {
    contact_frames.push_back(e.contactFrameId());
  }
  return contact_frames;
}


inline std::vector<std::string> Robot::contactFrameNames() const {
  std::vector<std::string> contact_frames;
  for (const auto& e : point_contacts_) {
    contact_frames.push_back(e.contactModelInfo().frame);
  }
  for (const auto& e : surface_contacts_) {
    contact_frames.push_back(e.contactModelInfo().frame);
  }
  return contact_frames;
}


inline std::vector<int> Robot::pointContactFrames() const {
  std::vector<int> contact_frames;
  for (const auto& e : point_contacts_) {
    contact_frames.push_back(e.contactFrameId());
  }
  return contact_frames;
}


inline std::vector<std::string> Robot::pointContactFrameNames() const {
  std::vector<std::string> contact_frames;
  for (const auto& e : point_contacts_) {
    contact_frames.push_back(e.contactModelInfo().frame);
  }
  return contact_frames;
}


inline std::vector<int> Robot::surfaceContactFrames() const {
  std::vector<int> contact_frames;
  for (const auto& e : surface_contacts_) {
    contact_frames.push_back(e.contactFrameId());
  }
  return contact_frames;
}


inline std::vector<std::string> Robot::surfaceContactFrameNames() const {
  std::vector<std::string> contact_frames;
  for (const auto& e : surface_contacts_) {
    contact_frames.push_back(e.contactModelInfo().frame);
  }
  return contact_frames;
}


inline ContactStatus Robot::createContactStatus() const {
  return ContactStatus(contactTypes(), contactFrameNames());
}


inline ImpulseStatus Robot::createImpulseStatus() const {
  return ImpulseStatus(contactTypes(), contactFrameNames());
}


inline const RobotProperties& Robot::robotProperties() const {
  return properties_;
}


inline void Robot::setRobotProperties(const RobotProperties& properties) {
  properties_ = properties;
  properties_.has_generalized_momentum_bias = false;
  if (properties_.generalized_momentum_bias.size() == dimv_) {
    properties_.has_generalized_momentum_bias 
        = !(properties_.generalized_momentum_bias.isZero());
  } 
}

} // namespace robotoc

#endif // ROBOTOC_ROBOT_HXX_ 