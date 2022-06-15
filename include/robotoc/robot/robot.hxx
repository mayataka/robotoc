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
  if (has_floating_base_) {
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
  int dimf_stack = 0;
  int point_contact_index = 0;
  int surface_contact_index = 0;
  for (int i=0; i<contact_frames_.size(); ++i) {
    switch (contact_types_[i]) {
      case ContactType::PointContact:
        if (contact_status.isContactActive(i)) {
          point_contacts_[point_contact_index].computeBaumgarteResidual(
              model_, data_, contact_status.contactPosition(i),
              (const_cast<Eigen::MatrixBase<VectorType>&>(baumgarte_residual))
                  .template segment<3>(dimf_stack));
          dimf_stack += 3;
        }
        point_contact_index += 1;
        break;
      case ContactType::SurfaceContact:
        if (contact_status.isContactActive(i)) {
          surface_contacts_[surface_contact_index].computeBaumgarteResidual(
              model_, data_, contact_status.contactPlacement(i),
              (const_cast<Eigen::MatrixBase<VectorType>&>(baumgarte_residual))
                  .template segment<6>(dimf_stack));
          dimf_stack += 6;
        }
        surface_contact_index += 1;
        break;
      default:
        break;
    }
  }
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
  int dimf_stack = 0;
  int point_contact_index = 0;
  int surface_contact_index = 0;
  for (int i=0; i<contact_frames_.size(); ++i) {
    switch (contact_types_[i]) {
      case ContactType::PointContact:
        if (contact_status.isContactActive(i)) {
          point_contacts_[point_contact_index].computeBaumgarteDerivatives(
              model_, data_, 
              (const_cast<Eigen::MatrixBase<MatrixType1>&>(baumgarte_partial_dq))
                  .block(dimf_stack, 0, 3, dimv_),
              (const_cast<Eigen::MatrixBase<MatrixType2>&>(baumgarte_partial_dv))
                  .block(dimf_stack, 0, 3, dimv_),
              (const_cast<Eigen::MatrixBase<MatrixType3>&>(baumgarte_partial_da))
                  .block(dimf_stack, 0, 3, dimv_));
          dimf_stack += 3;
        }
        point_contact_index += 1;
        break;
      case ContactType::SurfaceContact: 
        if (contact_status.isContactActive(i)) {
          surface_contacts_[surface_contact_index].computeBaumgarteDerivatives(
              model_, data_, 
              (const_cast<Eigen::MatrixBase<MatrixType1>&>(baumgarte_partial_dq))
                  .block(dimf_stack, 0, 6, dimv_),
              (const_cast<Eigen::MatrixBase<MatrixType2>&>(baumgarte_partial_dv))
                  .block(dimf_stack, 0, 6, dimv_),
              (const_cast<Eigen::MatrixBase<MatrixType3>&>(baumgarte_partial_da))
                  .block(dimf_stack, 0, 6, dimv_));
          dimf_stack += 6;
        }
        surface_contact_index += 1;
        break;
      default:
        break;
    }
  }
}


template <typename VectorType>
inline void Robot::computeImpulseVelocityResidual(
    const ImpulseStatus& impulse_status, 
    const Eigen::MatrixBase<VectorType>& velocity_residual) const {
  assert(velocity_residual.size() == impulse_status.dimi());
  int dimf_stack = 0;
  int point_contact_index = 0;
  int surface_contact_index = 0;
  for (int i=0; i<contact_frames_.size(); ++i) {
    switch (contact_types_[i]) {
      case ContactType::PointContact:
        if (impulse_status.isImpulseActive(i)) {
          point_contacts_[point_contact_index].computeContactVelocityResidual(
              model_, data_, 
              (const_cast<Eigen::MatrixBase<VectorType>&>(velocity_residual))
                  .template segment<3>(dimf_stack));
          dimf_stack += 3;
        }
        point_contact_index += 1;
        break;
      case ContactType::SurfaceContact:
        if (impulse_status.isImpulseActive(i)) {
          surface_contacts_[surface_contact_index].computeContactVelocityResidual(
              model_, data_, 
              (const_cast<Eigen::MatrixBase<VectorType>&>(velocity_residual))
                  .template segment<6>(dimf_stack));
          dimf_stack += 6;
        }
        surface_contact_index += 1;
        break;
      default:
        break;
    }
  }
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
  int dimf_stack = 0;
  int point_contact_index = 0;
  int surface_contact_index = 0;
  for (int i=0; i<contact_frames_.size(); ++i) {
    switch (contact_types_[i]) {
      case ContactType::PointContact:
        if (impulse_status.isImpulseActive(i)) {
          point_contacts_[point_contact_index].computeContactVelocityDerivatives(
              model_, data_, 
              (const_cast<Eigen::MatrixBase<MatrixType1>&>(velocity_partial_dq))
                  .block(dimf_stack, 0, 3, dimv_),
              (const_cast<Eigen::MatrixBase<MatrixType2>&>(velocity_partial_dv))
                  .block(dimf_stack, 0, 3, dimv_));
          dimf_stack += 3;
        }
        point_contact_index += 1;
        break;
      case ContactType::SurfaceContact:
        if (impulse_status.isImpulseActive(i)) {
          surface_contacts_[surface_contact_index].computeContactVelocityDerivatives(
              model_, data_, 
              (const_cast<Eigen::MatrixBase<MatrixType1>&>(velocity_partial_dq))
                  .block(dimf_stack, 0, 6, dimv_),
              (const_cast<Eigen::MatrixBase<MatrixType2>&>(velocity_partial_dv))
                  .block(dimf_stack, 0, 6, dimv_));
          dimf_stack += 6;
        }
        surface_contact_index += 1;
        break;
      default:
        break;
    }
  }
}


template <typename VectorType>
inline void Robot::computeContactPositionResidual(
    const ImpulseStatus& impulse_status, 
    const Eigen::MatrixBase<VectorType>& position_residual) {
  assert(position_residual.size() == impulse_status.dimi());
  int dimf_stack = 0;
  int point_contact_index = 0;
  int surface_contact_index = 0;
  for (int i=0; i<contact_frames_.size(); ++i) {
    switch (contact_types_[i]) {
      case ContactType::PointContact:
        if (impulse_status.isImpulseActive(i)) {
          point_contacts_[point_contact_index].computeContactPositionResidual(
              model_, data_, impulse_status.contactPosition(i),
              (const_cast<Eigen::MatrixBase<VectorType>&>(position_residual))
                  .template segment<3>(dimf_stack));
          dimf_stack += 3;
        }
        point_contact_index += 1;
        break;
      case ContactType::SurfaceContact:
        if (impulse_status.isImpulseActive(i)) {
          surface_contacts_[surface_contact_index].computeContactPositionResidual(
              model_, data_, impulse_status.contactPlacement(i),
              (const_cast<Eigen::MatrixBase<VectorType>&>(position_residual))
                  .template segment<6>(dimf_stack));
          dimf_stack += 6;
        }
        surface_contact_index += 1;
        break;
      default:
        break;
    }
  }
}


template <typename MatrixType>
inline void Robot::computeContactPositionDerivative(
    const ImpulseStatus& impulse_status, 
    const Eigen::MatrixBase<MatrixType>& position_partial_dq) {
  assert(position_partial_dq.rows() == impulse_status.dimi());
  assert(position_partial_dq.cols() == dimv_);
  int dimf_stack = 0;
  int point_contact_index = 0;
  int surface_contact_index = 0;
  for (int i=0; i<contact_frames_.size(); ++i) {
    switch (contact_types_[i]) {
      case ContactType::PointContact:
        if (impulse_status.isImpulseActive(i)) {
          point_contacts_[point_contact_index].computeContactPositionDerivative(
              model_, data_, 
            (const_cast<Eigen::MatrixBase<MatrixType>&>(position_partial_dq))
                .block(dimf_stack, 0, 3, dimv_));
          dimf_stack += 3;
        }
        point_contact_index += 1;
        break;
    case ContactType::SurfaceContact:
      if (impulse_status.isImpulseActive(i)) {
        surface_contacts_[surface_contact_index].computeContactPositionDerivative(
            model_, data_, 
          (const_cast<Eigen::MatrixBase<MatrixType>&>(position_partial_dq))
              .block(dimf_stack, 0, 6, dimv_));
        dimf_stack += 6;
      }
      surface_contact_index += 1;
      break;
    default:
      break;
    }
  }
}


inline void Robot::setContactForces(const ContactStatus& contact_status, 
                                    const std::vector<Vector6d>& f) {
  assert(f.size() == max_num_contacts_);
  int point_contact_index = 0;
  int surface_contact_index = 0;
  for (int i=0; i<contact_types_.size(); ++i) {
    switch (contact_types_[i]) {
      case ContactType::PointContact:
        if (contact_status.isContactActive(i)) {
          point_contacts_[point_contact_index].computeJointForceFromContactForce(
              f[i].template head<3>(), fjoint_);
        }
        else {
          point_contacts_[point_contact_index].computeJointForceFromContactForce(
              Eigen::Vector3d::Zero(), fjoint_);
        }
        ++point_contact_index;
        break;
      case ContactType::SurfaceContact:
        if (contact_status.isContactActive(i)) {
          surface_contacts_[surface_contact_index].computeJointForceFromContactWrench(
              f[i], fjoint_);
        }
        else {
          surface_contacts_[surface_contact_index].computeJointForceFromContactWrench(
              Vector6d::Zero(), fjoint_);
        }
        ++surface_contact_index;
        break;
      default:
        break;
    }
  }
}


inline void Robot::setImpulseForces(const ImpulseStatus& impulse_status, 
                                    const std::vector<Vector6d>& f) {
  assert(f.size() == max_num_contacts_);
  int point_contact_index = 0;
  int surface_contact_index = 0;
  for (int i=0; i<contact_types_.size(); ++i) {
    switch (contact_types_[i]) {
      case ContactType::PointContact:
        if (impulse_status.isImpulseActive(i)) {
          point_contacts_[point_contact_index].computeJointForceFromContactForce(
              f[i].template head<3>(), fjoint_);
        }
        else {
          point_contacts_[point_contact_index].computeJointForceFromContactForce(
              Eigen::Vector3d::Zero(), fjoint_);
        }
        ++point_contact_index;
        break;
      case ContactType::SurfaceContact:
        if (impulse_status.isImpulseActive(i)) {
          surface_contacts_[surface_contact_index].computeJointForceFromContactWrench(
              f[i], fjoint_);
        }
        else {
          surface_contacts_[surface_contact_index].computeJointForceFromContactWrench(
              Vector6d::Zero(), fjoint_);
        }
        ++surface_contact_index;
        break;
      default:
        break;
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
  if (contact_inv_damping_ > 0.) {
    data_.JMinvJt.diagonal().array() += contact_inv_damping_;
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
  if (has_floating_base_) {
    q_min.head(7) = - Eigen::VectorXd::Ones(7);
    q_max.head(7) = Eigen::VectorXd::Ones(7);
  }
  q_min.tail(dimu_) = lower_joint_position_limit_;
  q_max.tail(dimu_) = upper_joint_position_limit_;
  return pinocchio::randomConfiguration(model_, q_min, q_max);
}


template <typename ConfigVectorType>
inline void Robot::normalizeConfiguration(
    const Eigen::MatrixBase<ConfigVectorType>& q) const {
  assert(q.size() == dimq_);
  if (has_floating_base_) {
    if (q.template segment<4>(3).squaredNorm() 
          <= std::numeric_limits<double>::epsilon()) {
      (const_cast<Eigen::MatrixBase<ConfigVectorType>&> (q)).coeffRef(3) = 1;
    }
    pinocchio::normalize(model_, 
                         const_cast<Eigen::MatrixBase<ConfigVectorType>&>(q));
  }
}


inline int Robot::frameId(const std::string& frame_name) const {
  try {
    if (!model_.existFrame(frame_name)) {
      throw std::invalid_argument(
          "Invalid argument: frame " + frame_name + " does not exit!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
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
  return has_floating_base_;
}


inline int Robot::maxNumContacts() const {
  return max_num_contacts_;
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
  return contact_types_[contact_index];
}


inline const std::vector<ContactType>& Robot::contactTypes() const {
  return contact_types_;
}


inline std::vector<int> Robot::contactFrames() const {
  return contact_frames_;
}


inline std::vector<std::string> Robot::contactFrameNames() const {
  return contact_frame_names_;
}


inline std::vector<int> Robot::pointContactFrames() const {
  std::vector<int> point_contact_frames;
  for (const auto& e : point_contacts_) {
    point_contact_frames.push_back(e.contact_frame_id());
  }
  return point_contact_frames;
}


inline std::vector<std::string> Robot::pointContactFrameNames() const {
  std::vector<std::string> point_contact_frame_names;
  for (const auto& e : point_contacts_) {
    point_contact_frame_names.push_back(frameName(e.contact_frame_id()));
  }
  return point_contact_frame_names;
}


inline std::vector<int> Robot::surfaceContactFrames() const {
  std::vector<int> surface_contact_frames;
  for (const auto& e : surface_contacts_) {
    surface_contact_frames.push_back(e.contact_frame_id());
  }
  return surface_contact_frames;
}


inline std::vector<std::string> Robot::surfaceContactFrameNames() const {
  std::vector<std::string> surface_contact_frame_names;
  for (const auto& e : surface_contacts_) {
    surface_contact_frame_names.push_back(frameName(e.contact_frame_id()));
  }
  return surface_contact_frame_names;
}


inline ContactStatus Robot::createContactStatus() const {
  return ContactStatus(contact_types_, contact_frame_names_);
}


inline ImpulseStatus Robot::createImpulseStatus() const {
  return ImpulseStatus(contact_types_, contact_frame_names_);
}

} // namespace robotoc

#endif // ROBOTOC_ROBOT_HXX_ 