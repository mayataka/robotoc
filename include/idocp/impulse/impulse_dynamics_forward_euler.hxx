#ifndef IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HXX_
#define IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HXX_

#include "idocp/impulse/impulse_dynamics_forward_euler.hpp"

#include <assert.h>

namespace idocp {

inline ImpulseDynamicsForwardEuler::ImpulseDynamicsForwardEuler(
    const Robot& robot) 
  : schur_complement_(robot.dimv(), robot.max_dimf()),
    dImD_dq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dImD_ddv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dImD_df_full_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    MJTJnv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                        robot.dimv()+robot.max_dimf())), 
    MJTJinvImDCqv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                              2*robot.dimv())), 
    ldvf_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    ImDC_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    dimv_(robot.dimv()),
    dimf_(0) {
}


inline ImpulseDynamicsForwardEuler::ImpulseDynamicsForwardEuler() 
  : schur_complement_(),
    dImD_dq_(),
    dImD_ddv_(),
    dImD_df_full_(),
    MJTJinv_full_(), 
    MJTJinvImDCqv_full_(), 
    ldvf_full_(),
    ImDC_full_(),
    dimv_(0),
    dimf_(0) {
}


inline ImpulseDynamicsForwardEuler::~ImpulseDynamicsForwardEuler() {
}


inline void ImpulseDynamicsForwardEuler::linearizeImpulseDynamics(
    Robot& robot, const ContactStatus& contact_status,  
    const ImpulseSplitSolution& s, ImpulseKKTMatrix& kkt_matrix, 
    ImpulseKKTResidual& kkt_residual) { 
  assert(contact_status.hasActiveContacts());
  setContactStatus(contact_status);
  linearizeInverseImpulseDynamics(robot, contact_status, s, kkt_residual);
  linearizeContactVelocityConstraint(robot, contact_status, kkt_matrix, 
                                     kkt_residual);
  // augment inverse dynamics constraint
  kkt_residual.lq().noalias() += dImD_dq_.transpose() * s.beta;
  kkt_residual.ldv.noalias() += dImD_ddv_.transpose() * s.beta;
  kkt_residual.lf().noalias() += dImD_df_().transpose() * s.beta;
  // augment contact constraint
  kkt_residual.lq().noalias()
      += kkt_matrix.Cq_contact_velocity().transpose() * s.mu_contact_velocity();
  kkt_residual.lv().noalias()
      += kkt_matrix.Cv_contact_velocity().transpose() * s.mu_contact_velocity();
  kkt_residual.ldv.noalias() 
      += kkt_matrix.Cv_contact_velocity().transpose() * s.mu_contact_velocity();
}


inline void ImpulseDynamicsForwardEuler::condenseImpulseDynamics(
    Robot& robot, const ContactStatus& contact_status, 
    const ImpulseSplitSolution& s, ImpulseKKTMatrix& kkt_matrix, 
    ImpulseKKTResidual& kkt_residual) {
  assert(contact_status.hasActiveContacts());
  setContactStatus(contact_status);
  linearizeInverseImpulseDynamics(robot, contact_status, s, kkt_residual);
  linearizeContactVelocityConstraint(robot, contact_status, kkt_matrix, 
                                     kkt_residual);
  schur_complement_.invertWithZeroBottomRightCorner(
      dimv_, dimf_, dImD_ddv_, kkt_matrix.Cv_contact_velocity().transpose(), 
      MJTJinv_());
  MJTJinvImDCqv_().topLeftCorner(dimv_, dimv_).noalias() = MJTJinv_().topLeftCorner(dimv_, dimv_) * dImD_dq_;
  MJTJinvImDCqv_().topLeftCorner(dimv_, dimv_).noalias() += MJTJinv_().topRightCorner(dimv_, dimf_) * kkt_matrix.Cq_contact_velocity();
  MJTJinvImDCqv_().topRightCorner(dimv_, dimv_).noalias() = MJTJinv_().topRightCorner(dimv_, dimf_) * kkt_matrix.Cv_contact_velocity();
  MJTJinvImDCqv_().bottomLeftCorner(dimf_, dimv_).noalias() = MJTJinv_().bottomLeftCorner(dimf_, dimv_) * dImd_dq_;
  MJTJinvImDCqv_().bottomLeftCorner(dimf_, dimv_).noalias() += MJTJinv_().bottomRightCorner(dimf_, dimf_) * kkt_matrix.Cq_contact_velocity();
  MJTJinvImDCqv_().bottomLeftCorner(dimf_, dimv_).noalias() = MJTJinv_().bottomRightCorner(dimf_, dimf_) * kkt_matrix.Cv_contact_velocity();
  ldvf_().head(dimv_) = kkt_residual.ldv;
  ldvf_().tail(dimf_) = kkt_residual.lf();
  ImDC_().head(dimv_) = kkt_residual.dv_res;
  ImDC_().tail(dimf_) = kkt_residual.C_contact_velocity();
  kkt_matrix.Qxx().noalias() 
      += MJTJinvCqv_().transpose() * ldvf_().asDiagonal() * MJTJinvCqv_();
  kkt_residual.lx().noalias()
      += MJTJinvCqv_().transpose() * ldvf_() * MJTJinv_() * ImDC_();
  kkt_matrix.Fqv = - MJTJinvCqv_().topLeftCorner(dimv_, dimv_);
  kkt_matrix.Fvv.setIdentity(dimv_, dimv_);
  kkt_matrix.Fvv.noalias() -= MJTJinvCqv_().topRightCorner(dimv_, dimv_);
}


inline void ImpulseDynamicsForwardEuler::computeCondensedDirection(
    const ImpulseKKTMatrix& kkt_matrix, const ImpulseKKTResidual& kkt_residual, 
    const ImpulseSplitDirection& d_next, ImpulseSplitDirection& d) {
  d.ddv.noalias() = - MJTJinvImDCqv_().topRows(dimv_) * d.dx() 
                    - MJTJinv_().topRows(dimv_) * ImDC_();
  d.df().noalias() = - MJTJinvImDCqv_().bottomRows(dimf_) * d.dx()
                     - MJTJinv_().bottomRows(dimf_) * ImDC_();
  d.beta.noalias() = - MJTJinvImDCqv_().bottomRows(dimf_) * d.dx()
                     - MJTJinv_().bottomRows(dimf_) * ImDC_();
  d.mu().noalias() = - MJTJinvImDCqv_().bottomRows(dimf_) * d.dx()
                     - MJTJinv_().bottomRows(dimf_) * ImDC_();

}


inline void ImpulseDynamicsForwardEuler::computeImpulseDynamicsResidual(
    Robot& robot, const ContactStatus& contact_status,  
    const ImpulseSplitSolution& s, ImpulseKKTResidual& kkt_residual) {
  robot.setContactForces(contact_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, kkt_residual.dv_res);
  robot.computeContactVelocityResidual(contact_status, 
                                       kkt_residual.C_contact_velocity());
}


inline void ImpulseDynamicsForwardEuler::linearizeInverseImpulseDynamics(
    Robot& robot, const ContactStatus& contact_status, 
    const ImpulseSplitSolution& s, ImpulseKKTResidual& kkt_residual) {
  robot.setContactForces(contact_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, kkt_residual.dv_res);
  robot.RNEAImpulseDerivatives(s.q, s.dv, dImD_dq_, dImD_ddv_);
  robot.dRNEAPartialdFext(contact_status, dImD_df_full_);
}


inline void ImpulseDynamicsForwardEuler::linearizeContactVelocityConstraint(
    Robot& robot, const ContactStatus& contact_status, 
    ImpulseKKTMatrix& kkt_matrix, ImpulseKKTResidual& kkt_residual) {
  robot.computeContactVelocityResidual(contact_status, 
                                       kkt_residual.C_contact_velocity());
  robot.computeContactVelocityDerivatives(contact_status, 
                                          kkt_matrix.Cq_contact_velocity(),
                                          kkt_matrix.Cv_contact_velocity());
}


inline double ImpulseDynamicsForwardEuler::l1NormImpulseDynamicsResidual(
    const ImpulseKKTResidual& kkt_residual) const {
  return (kkt_residual.u_res.lpNorm<1>() 
          + kkt_residual.C().lpNorm<1>());
}


inline double ImpulseDynamicsForwardEuler::squaredNormImpulseDynamicsResidual(
    const ImpulseKKTResidual& kkt_residual) const {
  return (kkt_residual.dv_res.squaredNorm() 
          + kkt_residual.C().squaredNorm());
}


inline void ImpulseDynamicsForwardEuler::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
}


inline Eigen::Block<Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic, true> 
ImpulseDynamicsForwardEuler::dImD_df_() {
  return dImD_df_full_.leftCols(dimf_);
}


inline Eigen::Block<Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic, true> 
ImpulseDynamicsForwardEuler::MJTJinv_() {
  return MJTJinv_full_.topLeftCorner(dimv_+dimf_, dimv_+dimf_);
}


inline Eigen::Block<Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic, true> 
ImpulseDynamicsForwardEuler::MJTJinvImDCqv_() {
  return MJTJinvImDCqv_full_.topLeftCorner(dimv_+dimf_, 2*dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsForwardEuler::ldvf_() {
  return ldvf_full_.head(dimv_+dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsForwardEuler::ImDC_() {
  return ImDC_full_.head(dimv_+dimf_);
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HXX_ 