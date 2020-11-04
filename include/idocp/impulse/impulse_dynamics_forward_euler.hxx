#ifndef IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HXX_
#define IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HXX_

#include "idocp/impulse/impulse_dynamics_forward_euler.hpp"

#include <assert.h>

namespace idocp {

inline ImpulseDynamicsForwardEuler::ImpulseDynamicsForwardEuler(
    const Robot& robot) 
  : data_(robot),
    dimv_(robot.dimv()),
    dimf_(0) {
}


inline ImpulseDynamicsForwardEuler::ImpulseDynamicsForwardEuler() 
  : data_(robot),
    dimv_(robot.dimv()),
    dimf_(0) {
}


inline ImpulseDynamicsForwardEuler::~ImpulseDynamicsForwardEuler() {
}


inline void ImpulseDynamicsForwardEuler::linearizeImpulseDynamics(
    Robot& robot, const ImpulseStatus& impulse_status,  
    const ImpulseSplitSolution& s, ImpulseKKTMatrix& kkt_matrix, 
    ImpulseKKTResidual& kkt_residual) { 
  assert(contact_status.hasActiveContacts());
  setImpulseStatus(impulse_status);
  linearizeInverseImpulseDynamics(robot, impulse_status, s, data_);
  linearizeImpulseVelocityConstraint(robot, impulse_status, data_);
  linearizeImpulsePositionConstraint(robot, impulse_status, kkt_matrix, 
                                     kkt_residual);
  // augment inverse impulse dynamics constraint
  kkt_residual.lq().noalias() += data_.dImDdq.transpose() * s.beta;
  kkt_residual.ldv.noalias() += data_.dImDddv.transpose() * s.beta;
  kkt_residual.lf().noalias() -= kkt_matrix.Cv() * s.beta;
  // augment impulse velocity constraint
  kkt_residual.lq().noalias() += data_.dCdq().transpose() * s.mu_stack();
  kkt_residual.lv().noalias() += data_.dCdv().transpose() * s.mu_stack();
  kkt_residual.ldv.noalias() += data_.dCdv().transpose() * s.mu_stack();
  // augment impulse position constraint
  kkt_residual.lq().noalias() += kkt_matrix.Pq().transpose() * s.xi_stack();
}


inline void ImpulseDynamicsForwardEuler::linearizeInverseImpulseDynamics(
    Robot& robot, const ImpulseStatus& impulse_status, 
    const ImpulseSplitSolution& s, ImpulseDynamicsData& data) {
  robot.setImpulseForces(impulse_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, data.ID());
  robot.RNEAImpulseDerivatives(s.q, s.dv, data.dImDdq(), data.dImDda);
}


inline void ImpulseDynamicsForwardEuler::linearizeImpulseVelocityConstraint(
    Robot& robot, const ImpulseStatus& impulse_status, 
    ImpulseDynamicsData& data) {
  robot.computeImpulseVelocityResidual(impulse_status, data.C());
  robot.computeImpulseVelocityDerivatives(impulse_status, data.dCdq(), 
                                          data.dCdv());
}


inline void ImpulseDynamicsForwardEuler::linearizeImpulsePositionConstraint(
      Robot& robot, const ImpulseStatus& impulse_status, 
      ImpulseKKTMatrix& kkt_matrix, ImpulseKKTResidual& kkt_residual) {
  robot.computeImpulsePositionResidual(impulse_status, kkt_residual.P());
  robot.computeImpulsePositionDerivative(impulse_status, kkt_matrix.Pq());
}


inline void ImpulseDynamicsForwardEuler::condenseImpulseDynamics(
    Robot& robot, const ImpulseStatus& impulse_status, 
    ImpulseKKTMatrix& kkt_matrix, ImpulseKKTResidual& kkt_residual) {
  robot.computeMJtJinv(data_.dImDddv, data_.dCdv(), data_.MJtJinv());
  condensing(robot, data_, kkt_matrix, kkt_residual);
}


inline void ImpulseDynamicsForwardEuler::condensing(
    Robot& robot, const ImpulseStatus& impulse_status, 
    ImpulseKKTMatrix& kkt_matrix, ImpulseKKTResidual& kkt_residual) {
  const int dimv = robot.dimv();
  data.MJtJinv_dImDCdqv().topLeftCorner(dimv_, dimv_).noalias() 
      = MJtJinv_().topLeftCorner(dimv_, dimv_) * dImD_dq_;
  MJtJinv_dImDCdqv_().topLeftCorner(dimv_, dimv_).noalias() 
      += MJtJinv_().topRightCorner(dimv_, dimf_) * kkt_matrix.Cq();
  MJtJinv_dImDCdqv_().topRightCorner(dimv_, dimv_).noalias() 
      = MJtJinv_().topRightCorner(dimv_, dimf_) * kkt_matrix.Cv();
  MJtJinv_dImDCdqv_().bottomLeftCorner(dimf_, dimv_).noalias() 
      = MJtJinv_().bottomLeftCorner(dimf_, dimv_) * dImD_dq_;
  MJtJinv_dImDCdqv_().bottomLeftCorner(dimf_, dimv_).noalias() 
      += MJtJinv_().bottomRightCorner(dimf_, dimf_) * kkt_matrix.Cq();
  MJtJinv_dImDCdqv_().bottomRightCorner(dimf_, dimv_).noalias() 
      = MJtJinv_().bottomRightCorner(dimf_, dimf_) * kkt_matrix.Cv();
  MJtJinv_ImDC_().noalias() = MJtJinv_().leftCols(dimv_) * kkt_residual.dv_res;
  MJtJinv_ImDC_().noalias() += MJtJinv_().rightCols(dimf_) * kkt_residual.C();
  Qdvq_condensed_().noalias() = (- kkt_matrix.Qdvdv.diagonal()).asDiagonal() 
                                 * MJtJinv_dImDCdqv_().topLeftCorner(dimv_, dimv_);
  Qdvv_condensed_().noalias() = (- kkt_matrix.Qdvdv.diagonal()).asDiagonal() 
                                 * MJtJinv_dImDCdqv_().topRightCorner(dimv_, dimv_);
  Qfq_condensed_().noalias() = (- kkt_matrix.Qff().diagonal()).asDiagonal() 
                                * MJtJinv_dImDCdqv_().bottomLeftCorner(dimf_, dimv_);
  Qfv_condensed_().noalias() = (- kkt_matrix.Qff().diagonal()).asDiagonal() 
                                * MJtJinv_dImDCdqv_().bottomRightCorner(dimf_, dimv_);
  ldv_condensed_() = kkt_residual.ldv;
  ldv_condensed_().array() -= kkt_matrix.Qdvdv.diagonal().array() 
                              * MJtJinv_ImDC_().head(dimv_).array();
  lf_condensed_() = - kkt_residual.lf();
  lf_condensed_().array() -= kkt_matrix.Qff().diagonal().array() 
                              * MJtJinv_ImDC_().tail(dimf_).array();
  kkt_matrix.Qxx().noalias()
      -= MJtJinv_dImDCdqv_().transpose() * Qdvfqv_condensed_();
  kkt_residual.lx().noalias()
      -= MJtJinv_dImDCdqv_().transpose() * ldvf_condensed_();
  kkt_matrix.Fvq = - MJtJinv_dImDCdqv_().topLeftCorner(dimv_, dimv_);
  kkt_matrix.Fvv = Eigen::MatrixXd::Identity(dimv_, dimv_) 
                    - MJtJinv_dImDCdqv_().topRightCorner(dimv_, dimv_);
  kkt_residual.Fv().noalias() -= MJtJinv_ImDC_().head(dimv_);
}


inline void ImpulseDynamicsForwardEuler::computeCondensedDirection(
    const ImpulseKKTMatrix& kkt_matrix, const ImpulseKKTResidual& kkt_residual, 
    const SplitDirection& d_next, ImpulseSplitDirection& d) {
  d.ddv.noalias() = - MJtJinv_dImDCdqv_().topRows(dimv_) * d.dx() 
                    - MJtJinv_ImDC_().head(dimv_);
  d.df().noalias() = MJtJinv_dImDCdqv_().bottomRows(dimf_) * d.dx()
                     + MJtJinv_ImDC_().tail(dimf_);
  ldvf_condensed_().noalias() += Qdvfqv_condensed_() * d.dx();
  ldv_condensed_().noalias() += d_next.dgmm();
  d.dbeta.noalias() = - MJtJinv_().topRows(dimv_) * ldvf_condensed_();
  d.dmu().noalias() = - MJtJinv_().bottomRows(dimf_) * ldvf_condensed_();
}


inline void ImpulseDynamicsForwardEuler::computeImpulseDynamicsResidual(
    Robot& robot, const ContactStatus& contact_status,  
    const ImpulseSplitSolution& s, ImpulseKKTResidual& kkt_residual) {
  robot.setContactForces(contact_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, kkt_residual.dv_res);
  robot.computeContactVelocityResidual(contact_status, kkt_residual.C());
}


inline void ImpulseDynamicsForwardEuler::linearizeInverseImpulseDynamics(
    Robot& robot, const ContactStatus& contact_status, 
    const ImpulseSplitSolution& s, ImpulseKKTResidual& kkt_residual) {
  robot.setContactForces(contact_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, kkt_residual.dv_res);
  robot.RNEAImpulseDerivatives(s.q, s.dv, dImD_dq_, dImD_ddv_);
}


inline void ImpulseDynamicsForwardEuler::linearizeContactConstraint(
    Robot& robot, const ContactStatus& contact_status, 
    ImpulseKKTMatrix& kkt_matrix, ImpulseKKTResidual& kkt_residual) {
  robot.computeContactVelocityResidual(contact_status, kkt_residual.C());
  robot.computeContactVelocityDerivatives(contact_status, kkt_matrix.Cq(),
                                          kkt_matrix.Cv());
}


inline double ImpulseDynamicsForwardEuler::l1NormImpulseDynamicsResidual(
    const ImpulseKKTResidual& kkt_residual) {
  return (kkt_residual.dv_res.lpNorm<1>() + kkt_residual.C().lpNorm<1>());
}


inline double ImpulseDynamicsForwardEuler::squaredNormImpulseDynamicsResidual(
    const ImpulseKKTResidual& kkt_residual) {
  return (kkt_residual.dv_res.squaredNorm() + kkt_residual.C().squaredNorm());
}


inline void ImpulseDynamicsForwardEuler::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseDynamicsForwardEuler::MJtJinv_() {
  return MJtJinv_full_.topLeftCorner(dimv_+dimf_, dimv_+dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsForwardEuler::MJtJinv_dImDCdqv_() {
  return MJtJinv_dImDCdqv_full_.topLeftCorner(dimv_+dimf_, 2*dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsForwardEuler::Qdvfqv_condensed_() {
  return Qdvfqv_condensed_full_.topLeftCorner(dimv_+dimf_, 2*dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsForwardEuler::Qdvq_condensed_() {
  return Qdvfqv_condensed_full_.topLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsForwardEuler::Qdvv_condensed_() {
  return Qdvfqv_condensed_full_.block(0, dimv_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsForwardEuler::Qfq_condensed_() {
  return Qdvfqv_condensed_full_.block(dimv_, 0, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsForwardEuler::Qfv_condensed_() {
  return Qdvfqv_condensed_full_.block(dimv_, dimv_, dimf_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsForwardEuler::MJtJinv_ImDC_() {
  return MJtJinv_ImDC_full_.head(dimv_+dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsForwardEuler::ldvf_condensed_() {
  return ldvf_condensed_full_.head(dimv_+dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsForwardEuler::ldv_condensed_() {
  return ldvf_condensed_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsForwardEuler::lf_condensed_() {
  return ldvf_condensed_full_.segment(dimv_, dimf_);
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HXX_ 