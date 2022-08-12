#include "robotoc/constraints/friction_cone.hpp"

#include <stdexcept>
#include <iostream>


namespace robotoc {

FrictionCone::FrictionCone(const Robot& robot, const double mu)
  : FrictionCone(robot, std::vector<double>(robot.maxNumContacts(), mu)) {
}


FrictionCone::FrictionCone(const Robot& robot, const std::vector<double>& mu)
  : ConstraintComponentBase(),
    dimv_(robot.dimv()),
    dimc_(5*robot.maxNumContacts()),
    max_num_contacts_(robot.maxNumContacts()),
    contact_frame_(robot.contactFrames()),
    contact_types_(robot.contactTypes()),
    mu_(mu),
    cone_(robot.maxNumContacts(), Eigen::MatrixXd::Zero(5, 3)) {
  if (mu.size() != robot.maxNumContacts()) {
    throw std::out_of_range(
        "[FrictionCone] invalid argument: mu.size() must be" + std::to_string(robot.maxNumContacts()) + "!");
  }
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    if (mu[i] <= 0) {
      throw std::out_of_range(
          "[FrictionCone] invalid argument: mu[" + std::to_string(i) + "] must be positive!");
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


FrictionCone::FrictionCone()
  : ConstraintComponentBase(),
    dimc_(0),
    max_num_contacts_(0),
    dimv_(0),
    contact_frame_(),
    contact_types_(),
    mu_(),
    cone_() {
}


FrictionCone::~FrictionCone() {
}


void FrictionCone::setFrictionCoefficient(const double mu) {
  setFrictionCoefficient(std::vector<double>(max_num_contacts_, mu));
}


void FrictionCone::setFrictionCoefficient(const std::vector<double>& mu) {
  if (mu.size() != max_num_contacts_) {
    throw std::out_of_range(
        "[FrictionCone] invalid argument: mu.size() must be" + std::to_string(max_num_contacts_) + "!");
  }
  for (int i=0; i<max_num_contacts_; ++i) {
    if (mu[i] <= 0) {
      throw std::out_of_range(
          "[FrictionCone] invalid argument: mu[" + std::to_string(i) + "] must be positive!");
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


bool FrictionCone::useKinematics() const {
  return true;
}


KinematicsLevel FrictionCone::kinematicsLevel() const {
  return KinematicsLevel::AccelerationLevel;
}


void FrictionCone::allocateExtraData(ConstraintComponentData& data) const {
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


bool FrictionCone::isFeasible(Robot& robot, const ContactStatus& contact_status,
                              ConstraintComponentData& data, 
                              const SplitSolution& s) const {
  robot.updateFrameKinematics(s.q);
  for (int i=0; i<max_num_contacts_; ++i) {
    if (contact_status.isContactActive(i)) {
      const int idx = 5*i;
      Eigen::VectorXd& fWi = fW(data, i);
      robot.transformFromLocalToWorld(contact_frame_[i], 
                                      s.f[i].template head<3>(), fWi);
      frictionConeResidual(mu_[i], fWi, contact_status.contactRotation(i), 
                           data.residual.template segment<5>(idx));
      if (data.residual.maxCoeff() > 0) {
        return false;
      }
    }
  }
  return true;
}


void FrictionCone::setSlack(Robot& robot, const ContactStatus& contact_status, 
                            ConstraintComponentData& data, 
                            const SplitSolution& s) const {
  robot.updateFrameKinematics(s.q);
  for (int i=0; i<max_num_contacts_; ++i) {
    const int idx = 5*i;
    Eigen::VectorXd& fWi = fW(data, i);
    robot.transformFromLocalToWorld(contact_frame_[i], 
                                    s.f[i].template head<3>(), fWi);
    frictionConeResidual(mu_[i], fWi, contact_status.contactRotation(i), 
                         data.residual.template segment<5>(idx));
    data.slack.template segment<5>(idx)
        = - data.residual.template segment<5>(idx);
  }
}


void FrictionCone::evalConstraint(Robot& robot, 
                                  const ContactStatus& contact_status, 
                                  ConstraintComponentData& data, 
                                  const SplitSolution& s) const {
  data.residual.setZero();
  data.cmpl.setZero();
  data.log_barrier = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    if (contact_status.isContactActive(i)) {
      const int idx = 5*i;
      // Contact force expressed in the world frame.
      Eigen::VectorXd& fWi = fW(data, i);
      robot.transformFromLocalToWorld(contact_frame_[i], 
                                      s.f[i].template head<3>(), fWi);
      frictionConeResidual(mu_[i], fWi, contact_status.contactRotation(i), 
                           data.residual.template segment<5>(idx));
      data.residual.template segment<5>(idx).noalias()
          += data.slack.template segment<5>(idx);
      computeComplementarySlackness<5>(data, idx);
      data.log_barrier += logBarrier(data.slack.template segment<5>(idx));
    }
  }
}


void FrictionCone::evalDerivatives(Robot& robot, 
                                   const ContactStatus& contact_status, 
                                   ConstraintComponentData& data, 
                                   const SplitSolution& s, 
                                   SplitKKTResidual& kkt_residual) const {
  int dimf_stack = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    if (contact_status.isContactActive(i)) {
      const int idx = 5*i;
      // Contact force expressed in the world frame.
      const Eigen::VectorXd& fWi = fW(data, i);
      // Friction cone in the local frame of the contact surface.
      Eigen::MatrixXd& cone_local_i = cone_local(data, i);
      cone_local_i.noalias() = cone_[i] * contact_status.contactRotation(i).transpose();
      // Jacobian of the contact force expressed in the world frame fWi 
      // with respect to the configuration q.
      Eigen::MatrixXd& dfWi_dq = dfW_dq(data, i);
      robot.getJacobianTransformFromLocalToWorld(contact_frame_[i], fWi, dfWi_dq);
      // Jacobian of the frition cone constraint with respect to the 
      // configuration, i.e., s.q.
      Eigen::MatrixXd& dgi_dq = dg_dq(data, i);
      dgi_dq.noalias() = cone_local_i * dfWi_dq.template topRows<3>();
      kkt_residual.lq().noalias()
          += dgi_dq.transpose() * data.dual.template segment<5>(idx);
      // Jacobian of the frition cone constraint with respect to the contact
      // force expressed in the local frame, i.e., s.f[i].
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


void FrictionCone::condenseSlackAndDual(const ContactStatus& contact_status, 
                                        ConstraintComponentData& data, 
                                        SplitKKTMatrix& kkt_matrix, 
                                        SplitKKTResidual& kkt_residual) const {
  data.cond.setZero();
  int dimf_stack = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    if (contact_status.isContactActive(i)) {
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


void FrictionCone::expandSlackAndDual(const ContactStatus& contact_status,
                                      ConstraintComponentData& data, 
                                      const SplitDirection& d) const {
  // Because data.slack(i) and data.dual(i) are always positive,  
  // positive data.dslack and data.ddual do not affect the step size 
  // determined by the fraction-to-boundary-rule.
  data.dslack.fill(1.0);
  data.ddual.fill(1.0);
  int dimf_stack = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    if (contact_status.isContactActive(i)) {
      const int idx = 5*i;
      Eigen::MatrixXd& dgi_dq = dg_dq(data, i);
      Eigen::MatrixXd& dgi_df = dg_df(data, i);
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


int FrictionCone::dimc() const {
  return dimc_;
}

} // namespace robotoc