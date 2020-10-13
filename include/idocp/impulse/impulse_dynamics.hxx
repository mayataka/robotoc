#ifndef IDOCP_IMPULSE_DYNAMICS_HXX_
#define IDOCP_IMPULSE_DYNAMICS_HXX_

#include "idocp/ocp/impulse_dynamics.hpp"

#include <assert.h>

namespace idocp {

inline ImpulseDynamics::ImpulseDynamics(const Robot& robot) 
  : lu_condensed_(Eigen::VectorXd::Zero(robot.dimv())),
    C_floating_base_(Eigen::VectorXd::Zero(robot.dim_passive())),
    dImD_dq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dImD_ddv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dImD_df_full_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    Quu_du_dq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Quu_du_dv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Quu_du_da_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Quu_du_df_full_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    has_floating_base_(robot.has_floating_base()),
    has_active_contacts_(false),
    dimf_(0) {
}


inline ImpulseDynamics::ImpulseDynamics() 
  : lu_condensed_(),
    C_floating_base_(),
    dImD_dq_(),
    dImD_ddv_(),
    dImD_df_full_(),
    Quu_du_dq_(),
    Quu_du_dv_(),
    Quu_du_da_(),
    Quu_du_df_full_(),
    has_floating_base_(false),
    has_active_contacts_(false),
    dimf_(0) {
}


inline ImpulseDynamics::~ImpulseDynamics() {
}


inline void ImpulseDynamics::linearizeImpulseDynamics(
    Robot& robot, const ContactStatus& contact_status,  
    const ImpulseSplitSolution& s, ImpulseKKTMatrix& kkt_matrix, 
    ImpulseKKTResidual& kkt_residual) { 
  assert(contact_status.hasActiveContacts());
  setContactStatus(contact_status);
  linearizeInverseImpulseDynamics(robot, contact_status, s, kkt_residual);
  // augment inverse dynamics constraint
  kkt_residual.la().noalias() += dImD_da_.transpose() * s.beta;
  kkt_residual.lf().noalias() += dImD_df_().transpose() * s.beta;
  kkt_residual.lq().noalias() += dImD_dq_.transpose() * s.beta;
  linearizeContactVelocityConstraint(robot, contact_status, kkt_matrix, 
                                     kkt_residual);
  // augment contact constraint
  kkt_residual.lq().noalias()
      += kkt_matrix.Cq().transpose() * s.mu_contacts();
  kkt_residual.lv().noalias()
      += kkt_matrix.Cv().transpose() * s.mu_contacts();
  kkt_residual.ldv.noalias() 
      += kkt_matrix.Cdv.transpose() * s.mu_contacts();
}


inline void ImpulseDynamics::condenseImpulseDynamics(
    Robot& robot, const ContactStatus& contact_status, 
    const ImpulseSplitSolution& s, ImpulseKKTMatrix& kkt_matrix, 
    ImpulseKKTResidual& kkt_residual) {
  assert(contact_status.hasActiveContacts());
  setContactStatus(contact_status);
  linearizeInverseImpulseDynamics(robot, contact_status, s, kkt_residual);
  lu_condensed_.noalias() 
      = kkt_residual.lu 
          + kkt_matrix.Quu.diagonal().asDiagonal() * kkt_residual.u_res;
  // condense Newton residual
  kkt_residual.la().noalias() = du_da_.transpose() * lu_condensed_;
  kkt_residual.lf().noalias() = du_df_().transpose() * lu_condensed_;
  kkt_residual.lq().noalias() = du_dq_.transpose() * lu_condensed_;
  kkt_residual.lv().noalias() = du_dv_.transpose() * lu_condensed_;
  kkt_residual.lu.noalias() -= dtau * s.beta;   
  // condense Hessian
  // Assume that Quu has only diagonal elements, hence we write as 
  // kkt_matrix.Quu.dianobal().asDiagonal() for efficiency.
  // Otherwise, simply use kkt_matrix.Quu instead.
  Quu_du_da_.noalias() = kkt_matrix.Quu.diagonal().asDiagonal() * du_da_;
  kkt_matrix.Qaa().noalias() = du_da_.transpose() * Quu_du_da_;
  if (has_active_contacts_) {
    Quu_du_df_().noalias() = kkt_matrix.Quu.diagonal().asDiagonal() * du_df_(); 
    kkt_matrix.Qaf().noalias() = du_da_.transpose() * Quu_du_df_(); 
  }
  Quu_du_dq_.noalias() = kkt_matrix.Quu.diagonal().asDiagonal() * du_dq_;
  Quu_du_dv_.noalias() = kkt_matrix.Quu.diagonal().asDiagonal() * du_dv_;
  kkt_matrix.Qaq().noalias() = du_da_.transpose() * Quu_du_dq_;
  kkt_matrix.Qav().noalias() = du_da_.transpose() * Quu_du_dv_;
  if (has_active_contacts_) {
    kkt_matrix.Qff().noalias() = du_df_().transpose() * Quu_du_df_();
    kkt_matrix.Qfq().noalias() = du_df_().transpose() * Quu_du_dq_;
    kkt_matrix.Qfv().noalias() = du_df_().transpose() * Quu_du_dv_;
  }
  kkt_matrix.Qqq().noalias() = du_dq_.transpose() * Quu_du_dq_;
  kkt_matrix.Qqv().noalias() = du_dq_.transpose() * Quu_du_dv_;
  kkt_matrix.Qvv().noalias() = du_dv_.transpose() * Quu_du_dv_;
  // augment contact constraint 
  linearizeContactConstraint(robot, contact_status, dtau, kkt_matrix, 
                              kkt_residual);
  kkt_residual.la().noalias() += kkt_matrix.Ca().transpose() * s.mu_stack();
  kkt_residual.lq().noalias() += kkt_matrix.Cq().transpose() * s.mu_stack();
  kkt_residual.lv().noalias() += kkt_matrix.Cv().transpose() * s.mu_stack();
  }
}


inline void ImpulseDynamics::computeCondensedDirection(
    const ImpulseKKTMatrix& kkt_matrix, const ImpulseKKTResidual& kkt_residual, 
    ImpulseSplitDirection& d) {
  assert(dtau > 0);
  d.du = kkt_residual.u_res;
  d.du.noalias() += du_dq_ * d.dq();
  d.du.noalias() += du_dv_ * d.dv();
  d.du.noalias() += du_da_ * d.da();
  if (has_active_contacts_) {
    d.du.noalias() += du_df_() * d.df();
  }
  d.dbeta.noalias() = (kkt_residual.lu  + kkt_matrix.Quu * d.du) / dtau;
  if (has_floating_base_) {
    d.dbeta.template head<kDimFloatingBase>().noalias() 
        += d.dmu().template head<kDimFloatingBase>();
  }
}


inline void ImpulseDynamics::computeImpulseDynamicsResidual(
    Robot& robot, const ContactStatus& contact_status,  
    const ImpulseSplitSolution& s, ImpulseKKTResidual& kkt_residual) {
  robot.setContactForces(contact_status, s.f);
  robot.RNEAImpulse(s.q, s.a, kkt_residual.u_res);
  robot.computeContactVelocityResidual(contact_status, 
                                       kkt_residual.C_contacts());
}


inline void ImpulseDynamics::linearizeInverseImpulseDynamics(
    Robot& robot, const ContactStatus& contact_status, 
    const ImpulseSplitSolution& s, ImpulseKKTResidual& kkt_residual) {
  robot.setContactForces(contact_status, s.f);
  robot.RNEAImpulse(s.q, s.a, kkt_residual.u_res);
  robot.RNEAImpulseDerivatives(s.q, s.a, dImD_dq_, dImD_ddv_);
  robot.dRNEAPartialdFext(contact_status, dImD_df_full_);
}


inline void ImpulseDynamics::linearizeContactVelocityConstraint(
    Robot& robot, const ContactStatus& contact_status, 
    ImpulseKKTMatrix& kkt_matrix, ImpulseKKTResidual& kkt_residual) {
  robot.computeContactVelocityResidual(contact_status, 
                                       kkt_residual.C_contacts());
  robot.computeContactVelocityDerivatives(contact_status, dContactVel_dq_,
                                          dContactVel_dv_);
}


inline double ImpulseDynamics::l1NormImpulseDynamicsResidual(
    const ImpulseKKTResidual& kkt_residual) const {
  return (kkt_residual.u_res.lpNorm<1>() 
          + kkt_residual.C_contacts().lpNorm<1>());
}


inline double ImpulseDynamics::squaredNormImpulseDynamicsResidual(
    const ImpulseKKTResidual& kkt_residual) const {
  return (kkt_residual.u_res.squaredNorm() 
          + kkt_residual.C_contacts().squaredNorm());
}


inline void ImpulseDynamics::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
}


inline Eigen::Block<Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic, true> 
ImpulseDynamics::dImD_df_() {
  return dImD_df_full_.leftCols(dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic, true> 
ImpulseDynamics::dImD_df_() const {
  return dImD_df_full_.leftCols(dimf_);
}


inline Eigen::Block<Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic, true> 
ImpulseDynamics::Quu_du_df_() {
  return Quu_du_df_full_.leftCols(dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic, true> 
ImpulseDynamics::Quu_du_df_() const {
  return Quu_du_df_full_.leftCols(dimf_);
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_DYNAMICS_HXX_ 