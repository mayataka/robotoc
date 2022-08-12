#include "robotoc/constraints/impulse_friction_cone.hpp"

#include <stdexcept>
#include <iostream>


namespace robotoc {

ImpulseFrictionCone::ImpulseFrictionCone(const Robot& robot, const double mu)
  : ImpulseFrictionCone(robot, std::vector<double>(robot.maxNumContacts(), mu)) {
}


ImpulseFrictionCone::ImpulseFrictionCone(const Robot& robot, const std::vector<double>& mu)
  : ImpulseConstraintComponentBase(),
    dimv_(robot.dimv()),
    dimc_(5*robot.maxNumContacts()),
    max_num_contacts_(robot.maxNumContacts()),
    contact_frame_(robot.contactFrames()),
    contact_types_(robot.contactTypes()),
    mu_(mu),
    cone_(robot.maxNumContacts(), Eigen::MatrixXd::Zero(5, 3)) {
  if (mu.size() != robot.maxNumContacts()) {
    throw std::out_of_range(
        "Invalid argument: mu.size() must be" + std::to_string(robot.maxNumContacts()) + "!");
  }
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    if (mu[i] <= 0) {
      throw std::out_of_range(
          "Invalid argument: mu[" + std::to_string(i) + "] must be positive!");
    }
  }
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    cone_[i] <<  0,  0, -1, 
                 1,  0, -(mu[i]/std::sqrt(2)),
                -1,  0, -(mu[i]/std::sqrt(2)),
                 0,  1, -(mu[i]/std::sqrt(2)),
                 0, -1, -(mu[i]/std::sqrt(2));
  }
}


ImpulseFrictionCone::ImpulseFrictionCone()
  : ImpulseConstraintComponentBase(),
    dimc_(0),
    max_num_contacts_(0),
    dimv_(0),
    contact_frame_(),
    contact_types_(),
    mu_(),
    cone_() {
}


ImpulseFrictionCone::~ImpulseFrictionCone() {
}


void ImpulseFrictionCone::setFrictionCoefficient(const double mu) {
  setFrictionCoefficient(std::vector<double>(max_num_contacts_, mu));
}


void ImpulseFrictionCone::setFrictionCoefficient(const std::vector<double>& mu) {
  if (mu.size() != max_num_contacts_) {
    throw std::out_of_range(
        "Invalid argument: mu.size() must be" + std::to_string(max_num_contacts_) + "!");
  }
  for (int i=0; i<max_num_contacts_; ++i) {
    if (mu[i] <= 0) {
      throw std::out_of_range(
          "Invalid argument: mu[" + std::to_string(i) + "] must be positive!");
    }
  }
  for (int i=0; i<max_num_contacts_; ++i) {
    mu_[i]  = mu[i];
    cone_[i] <<  0,  0, -1, 
                 1,  0, -(mu[i]/std::sqrt(2)),
                -1,  0, -(mu[i]/std::sqrt(2)),
                 0,  1, -(mu[i]/std::sqrt(2)),
                 0, -1, -(mu[i]/std::sqrt(2));
  }
}


KinematicsLevel ImpulseFrictionCone::kinematicsLevel() const {
  return KinematicsLevel::AccelerationLevel;
}


void ImpulseFrictionCone::allocateExtraData(
    ConstraintComponentData& data) const {
  data.r.clear();
  for (int i=0; i<max_num_contacts_; ++i) {
    data.r.push_back(Eigen::VectorXd::Zero(3)); // fWi
  }
  for (int i=0; i<max_num_contacts_; ++i) {
    data.r.push_back(Eigen::VectorXd::Zero(5)); // ri
  }
  data.J.clear();
  for (int i=0; i<max_num_contacts_; ++i) {
    data.J.push_back(Eigen::MatrixXd::Zero(5, dimv_)); // dgi_dq
  }
  for (int i=0; i<max_num_contacts_; ++i) {
    data.J.push_back(Eigen::MatrixXd::Zero(5, 3)); // dgi_df
  }
  for (int i=0; i<max_num_contacts_; ++i) {
    data.J.push_back(Eigen::MatrixXd::Zero(6, dimv_)); // dfWi_dq
  }
  for (int i=0; i<max_num_contacts_; ++i) {
    data.J.push_back(Eigen::MatrixXd::Zero(5, 3)); // r_dgi_df
  }
  for (int i=0; i<max_num_contacts_; ++i) {
    data.J.push_back(Eigen::MatrixXd::Zero(5, 3)); // cone_local
  }
}


bool ImpulseFrictionCone::isFeasible(Robot& robot, 
                                     const ImpulseStatus& impulse_status,
                                     ConstraintComponentData& data, 
                                     const ImpulseSplitSolution& s) const {
  robot.updateFrameKinematics(s.q);
  for (int i=0; i<max_num_contacts_; ++i) {
    if (impulse_status.isImpulseActive(i)) {
      const int idx = 5*i;
      Eigen::VectorXd& fWi = fW(data, i);
      robot.transformFromLocalToWorld(contact_frame_[i], 
                                      s.f[i].template head<3>(), fWi);
      frictionConeResidual(mu_[i], fWi, impulse_status.contactRotation(i), 
                           data.residual.template segment<5>(idx));
      if (data.residual.maxCoeff() > 0) {
        return false;
      }
    }
  }
  return true;
}


void ImpulseFrictionCone::setSlack(Robot& robot, 
                                   const ImpulseStatus& impulse_status,
                                   ConstraintComponentData& data, 
                                   const ImpulseSplitSolution& s) const {
  robot.updateFrameKinematics(s.q);
  for (int i=0; i<max_num_contacts_; ++i) {
    const int idx = 5*i;
    Eigen::VectorXd& fWi = fW(data, i);
    robot.transformFromLocalToWorld(contact_frame_[i], 
                                    s.f[i].template head<3>(), fWi);
    frictionConeResidual(mu_[i], fWi, impulse_status.contactRotation(i), 
                         data.residual.template segment<5>(idx));
    data.slack.template segment<5>(idx)
        = - data.residual.template segment<5>(idx);
  }
}


void ImpulseFrictionCone::evalConstraint(Robot& robot, 
                                         const ImpulseStatus& impulse_status,
                                         ConstraintComponentData& data, 
                                         const ImpulseSplitSolution& s) const {
  data.residual.setZero();
  data.cmpl.setZero();
  data.log_barrier = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    if (impulse_status.isImpulseActive(i)) {
      const int idx = 5*i;
      // Contact force expressed in the world frame.
      Eigen::VectorXd& fWi = fW(data, i);
      robot.transformFromLocalToWorld(contact_frame_[i], 
                                      s.f[i].template head<3>(), fWi);
      frictionConeResidual(mu_[i], fWi, impulse_status.contactRotation(i), 
                           data.residual.template segment<5>(idx));
      data.residual.template segment<5>(idx).noalias()
          += data.slack.template segment<5>(idx);
      computeComplementarySlackness<5>(data, idx);
      data.log_barrier += logBarrier(data.slack.template segment<5>(idx));
    }
  }
}


void ImpulseFrictionCone::evalDerivatives(
    Robot& robot, const ImpulseStatus& impulse_status,
    ConstraintComponentData& data, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  int dimf_stack = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    if (impulse_status.isImpulseActive(i)) {
      const int idx = 5*i;
      // Contact force expressed in the world frame.
      const Eigen::VectorXd& fWi = fW(data, i);
      // Friction cone in the local frame of the contact surface.
      Eigen::MatrixXd& cone_local_i = cone_local(data, i);
      cone_local_i.noalias() = cone_[i] * impulse_status.contactRotation(i).transpose();
      // Jacobian of the contact force expressed in the world frame fWi 
      // with respect to the configuration q.
      Eigen::MatrixXd& dfWi_dq = dfW_dq(data, i);
      robot.getJacobianTransformFromLocalToWorld(contact_frame_[i], fWi, dfWi_dq);
      // Jacobian of the frition cone constraint with respect to the 
      // configuration q.
      Eigen::MatrixXd& dgi_dq = dg_dq(data, i);
      dgi_dq.noalias() = cone_local_i * dfWi_dq.template topRows<3>();
      kkt_residual.lq().noalias()
          += dgi_dq.transpose() * data.dual.template segment<5>(idx);
      // Jacobian of the frition cone constraint with respect to the contact
      // force expressed in the local frame.
      Eigen::MatrixXd& dgi_df = dg_df(data, i);
      dgi_df.noalias() = cone_local_i * robot.frameRotation(contact_frame_[i]);
      kkt_residual.lf().template segment<3>(dimf_stack).noalias()
          += dgi_df.transpose() * data.dual.template segment<5>(idx);
      switch (contact_types_[i]) {
        case ContactType::PointContact:
          dimf_stack += 3;
          break;
        case ContactType::SurfaceContact:
          dimf_stack += 6;
          break;
        default:
          break;
      }
    }
  }
}


void ImpulseFrictionCone::condenseSlackAndDual(
    const ImpulseStatus& impulse_status, ConstraintComponentData& data, 
    ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  data.cond.setZero();
  int dimf_stack = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    if (impulse_status.isImpulseActive(i)) {
      const int idx = 5*i;
      computeCondensingCoeffcient<5>(data, idx);
      const Vector5d& condi = data.cond.template segment<5>(idx);
      const Eigen::MatrixXd& dgi_dq = dg_dq(data, i);
      const Eigen::MatrixXd& dgi_df = dg_df(data, i);
      kkt_residual.lq().noalias() += dgi_dq.transpose() * condi;
      kkt_residual.lf().template segment<3>(dimf_stack).noalias()
          += dgi_df.transpose() * condi;
      Eigen::MatrixXd& dfWi_dq = dfW_dq(data, i);
      Eigen::MatrixXd& r_dgi_df = r_dg_df(data, i);
      Eigen::VectorXd& ri = r(data, i);
      ri.array() = data.dual.template segment<5>(idx).array() 
                    / data.slack.template segment<5>(idx).array();
      dfWi_dq.template topRows<5>().noalias() = ri.asDiagonal() * dgi_dq;
      r_dgi_df.noalias() = ri.asDiagonal() * dgi_df;
      kkt_matrix.Qqq().noalias()
          += dgi_dq.transpose() * dfWi_dq.template topRows<5>();
      kkt_matrix.Qqf().template middleCols<3>(dimf_stack).noalias()
          += dgi_dq.transpose() * r_dgi_df; 
      kkt_matrix.Qff().template block<3, 3>(dimf_stack, dimf_stack).noalias()
          += dgi_df.transpose() * r_dgi_df;
      switch (contact_types_[i]) {
        case ContactType::PointContact:
          dimf_stack += 3;
          break;
        case ContactType::SurfaceContact:
          dimf_stack += 6;
          break;
        default:
          break;
      }
    }
  }
}


void ImpulseFrictionCone::expandSlackAndDual(
     const ImpulseStatus& impulse_status, ConstraintComponentData& data, 
    const ImpulseSplitDirection& d) const {
  // Because data.slack(i) and data.dual(i) are always positive,  
  // positive data.dslack and data.ddual do not affect the step size 
  // determined by the fraction-to-boundary-rule.
  data.dslack.fill(1.0);
  data.ddual.fill(1.0);
  int dimf_stack = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    if (impulse_status.isImpulseActive(i)) {
      const int idx = 5*i;
      const Eigen::MatrixXd& dgi_dq = dg_dq(data, i);
      const Eigen::MatrixXd& dgi_df = dg_df(data, i);
      data.dslack.template segment<5>(idx).noalias()
          = - dgi_dq * d.dq() - dgi_df * d.df().template segment<3>(dimf_stack) 
            - data.residual.template segment<5>(idx);
      computeDualDirection<5>(data, idx);
      switch (contact_types_[i]) {
        case ContactType::PointContact:
          dimf_stack += 3;
          break;
        case ContactType::SurfaceContact:
          dimf_stack += 6;
          break;
        default:
          break;
      }
    }
  }
}


int ImpulseFrictionCone::dimc() const {
  return dimc_;
}

} // namespace robotoc