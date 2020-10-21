#ifndef IDOCP_IMPULSE_DYNAMICS_BACKWARD_EULER_HXX
#define IDOCP_IMPULSE_DYNAMICS_BACKWARD_EULER_HXX

#include "idocp/impulse/impulse_dynamics_backward_euler.hpp"

#include "pinocchio/algorithm/cholesky.hpp"

#include <assert.h>

namespace idocp {

inline ImpulseDynamicsBackwardEuler::ImpulseDynamicsBackwardEuler(
    const Robot& robot) 
  : dImD_dq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dImD_ddv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dImD_df_full_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    Minv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    MinvImDq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    MinvImDf_full_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    Qdvq_condensed_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qdvf_condensed_full_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    MinvImD_(Eigen::VectorXd::Zero(robot.dimv())),
    ldv_condensed_(Eigen::VectorXd::Zero(robot.dimv())),
    dimv_(robot.dimv()),
    dimf_(0) {
}


inline ImpulseDynamicsBackwardEuler::ImpulseDynamicsBackwardEuler() 
  : dImD_dq_(), 
    dImD_ddv_(), 
    dImD_df_full_(), 
    Minv_(), 
    MinvImDq_(), 
    MinvImDf_full_(), 
    Qdvq_condensed_(), 
    Qdvf_condensed_full_(),
    MinvImD_(), 
    ldv_condensed_(),
    dimv_(0), 
    dimf_(0) {
}


inline ImpulseDynamicsBackwardEuler::~ImpulseDynamicsBackwardEuler() {
}


inline void ImpulseDynamicsBackwardEuler::linearizeImpulseDynamics(
    Robot& robot, const ContactStatus& contact_status,  
    const ImpulseSplitSolution& s, ImpulseKKTMatrix& kkt_matrix, 
    ImpulseKKTResidual& kkt_residual) { 
  assert(contact_status.hasActiveContacts());
  setContactStatus(contact_status);
  linearizeInverseImpulseDynamics(robot, contact_status, s, kkt_residual);
  linearizeContactConstraint(robot, contact_status, kkt_matrix, kkt_residual);
  // augment inverse dynamics constraint
  kkt_residual.lq().noalias() += dImD_dq_.transpose() * s.beta;
  kkt_residual.ldv.noalias() += dImD_ddv_.transpose() * s.beta;
  kkt_residual.lf().noalias() += dImD_df_().transpose() * s.beta;
  // augment contact constraint
  kkt_residual.lq().noalias() += kkt_matrix.Cq().transpose() * s.mu_stack();
  kkt_residual.lv().noalias() += kkt_matrix.Cv().transpose() * s.mu_stack();
  kkt_residual.ldv.noalias() += kkt_matrix.Cv().transpose() * s.mu_stack();
}


inline void ImpulseDynamicsBackwardEuler::condenseImpulseDynamics(
    Robot& robot, const ContactStatus& contact_status, 
    const ImpulseSplitSolution& s, ImpulseKKTMatrix& kkt_matrix, 
    ImpulseKKTResidual& kkt_residual) {
  assert(contact_status.hasActiveContacts());
  linearizeImpulseDynamics(robot, contact_status, s, kkt_matrix, kkt_residual);
  robot.computeMinv(dImD_ddv_, Minv_);
  MinvImDq_.noalias() = Minv_ * dImD_dq_;
  MinvImDf_().noalias() = Minv_ * dImD_df_();
  MinvImD_.noalias() = Minv_ * kkt_residual.dv_res;
  Qdvq_condensed_.noalias() 
      = (- kkt_matrix.Qdvdv.diagonal()).asDiagonal() * MinvImDq_;
  Qdvf_condensed_().noalias() 
      = (- kkt_matrix.Qdvdv.diagonal()).asDiagonal() * MinvImDf_();
  ldv_condensed_.noalias() 
      = kkt_residual.ldv - kkt_matrix.Qdvdv.diagonal().asDiagonal() * MinvImD_;
  kkt_matrix.Qqq().noalias() -= MinvImDq_.transpose() * Qdvq_condensed_;
  kkt_matrix.Qqf().noalias() -= MinvImDq_.transpose() * Qdvf_condensed_();
  kkt_residual.lq().noalias() -= MinvImDq_.transpose() * ldv_condensed_;
  kkt_matrix.Qff().noalias() -= MinvImDf_().transpose() * Qdvf_condensed_();
  kkt_residual.lf().noalias() -= MinvImDf_().transpose() * ldv_condensed_;
  kkt_matrix.Fvq = - MinvImDq_;
  kkt_matrix.Fvf() = - MinvImDf_();
  kkt_residual.Fv().noalias() -= MinvImD_;
}


inline void ImpulseDynamicsBackwardEuler::computeCondensedDirection(
    const ImpulseKKTMatrix& kkt_matrix, const ImpulseKKTResidual& kkt_residual, 
    ImpulseSplitDirection& d) const {
  d.ddv = - MinvImDq_ * d.dq() - MinvImDf_() * d.df() - MinvImD_;
  d.dbeta = - Minv_ * Qdvq_condensed_ * d.dq() 
            - Minv_ * Qdvf_condensed_() * d.df()
            - Minv_ * d.dgmm() - Minv_ * ldv_condensed_;
}


inline void ImpulseDynamicsBackwardEuler::computeImpulseDynamicsResidual(
    Robot& robot, const ContactStatus& contact_status,  
    const ImpulseSplitSolution& s, ImpulseKKTResidual& kkt_residual) {
  robot.setContactForces(contact_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, kkt_residual.dv_res);
  robot.computeContactVelocityResidual(contact_status, kkt_residual.C());
}


inline void ImpulseDynamicsBackwardEuler::linearizeInverseImpulseDynamics(
    Robot& robot, const ContactStatus& contact_status, 
    const ImpulseSplitSolution& s, ImpulseKKTResidual& kkt_residual) {
  robot.setContactForces(contact_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, kkt_residual.dv_res);
  robot.RNEAImpulseDerivatives(s.q, s.dv, dImD_dq_, dImD_ddv_);
  robot.dRNEAPartialdFext(contact_status, dImD_df_full_);
}


inline void ImpulseDynamicsBackwardEuler::linearizeContactConstraint(
    Robot& robot, const ContactStatus& contact_status, 
    ImpulseKKTMatrix& kkt_matrix, ImpulseKKTResidual& kkt_residual) {
  robot.computeContactVelocityResidual(contact_status, kkt_residual.C());
  robot.computeContactVelocityDerivatives(contact_status, kkt_matrix.Cq(),
                                          kkt_matrix.Cv());
}


inline double ImpulseDynamicsBackwardEuler::l1NormImpulseDynamicsResidual(
    const ImpulseKKTResidual& kkt_residual) {
  return (kkt_residual.dv_res.lpNorm<1>() + kkt_residual.C().lpNorm<1>());
}


inline double ImpulseDynamicsBackwardEuler::squaredNormImpulseDynamicsResidual(
    const ImpulseKKTResidual& kkt_residual) {
  return (kkt_residual.dv_res.squaredNorm() + kkt_residual.C().squaredNorm());
}


inline void ImpulseDynamicsBackwardEuler::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseDynamicsBackwardEuler::dImD_df_() {
  return dImD_df_full_.topLeftCorner(dimv_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsBackwardEuler::dImD_df_() const {
  return dImD_df_full_.topLeftCorner(dimv_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseDynamicsBackwardEuler::MinvImDf_() {
  return MinvImDf_full_.topLeftCorner(dimv_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsBackwardEuler::MinvImDf_() const {
  return MinvImDf_full_.topLeftCorner(dimv_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsBackwardEuler::Qdvf_condensed_() {
  return Qdvf_condensed_full_.topLeftCorner(dimv_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsBackwardEuler::Qdvf_condensed_() const {
  return Qdvf_condensed_full_.topLeftCorner(dimv_, dimf_);
}


} // namespace idocp 

#endif // IDOCP_IMPULSE_DYNAMICS_BACKWARD_EULER_HXX 