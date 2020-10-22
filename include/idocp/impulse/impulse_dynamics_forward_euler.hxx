#ifndef IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HXX_
#define IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HXX_

#include "idocp/impulse/impulse_dynamics_forward_euler.hpp"

#include <assert.h>

namespace idocp {

inline ImpulseDynamicsForwardEuler::ImpulseDynamicsForwardEuler(
    const Robot& robot) 
  : dImD_dq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dImD_ddv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dImD_df_full_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    MJtJinv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                        robot.dimv()+robot.max_dimf())), 
    MJtJinvImDCqv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                              2*robot.dimv())), 
    Qdvq_condensed_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qdvv_condensed_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qfq_condensed_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
    Qfv_condensed_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
    MJtJinvImDC_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    ldv_condensed_(Eigen::VectorXd::Zero(robot.dimv())),
    lf_condensed_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    dimv_(robot.dimv()),
    dimf_(0) {
}


inline ImpulseDynamicsForwardEuler::ImpulseDynamicsForwardEuler() 
  : dImD_dq_(),
    dImD_ddv_(),
    dImD_df_full_(),
    MJtJinv_full_(), 
    MJtJinvImDCqv_full_(), 
    Qdvq_condensed_(),
    Qdvv_condensed_(),
    Qfq_condensed_full_(),
    Qfv_condensed_full_(),
    MJtJinvImDC_full_(),
    ldv_condensed_(),
    lf_condensed_full_(),
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


inline void ImpulseDynamicsForwardEuler::condenseImpulseDynamics(
    Robot& robot, const ContactStatus& contact_status, 
    const ImpulseSplitSolution& s, ImpulseKKTMatrix& kkt_matrix, 
    ImpulseKKTResidual& kkt_residual) {
  assert(contact_status.hasActiveContacts());
  linearizeImpulseDynamics(robot, contact_status, s, kkt_matrix, kkt_residual);
  robot.computeMJtJinv(contact_status, dImD_ddv_, kkt_matrix.Cv(), MJtJinv_());
  MJtJinvImDCqv_().topLeftCorner(dimv_, dimv_).noalias() 
      = MJtJinv_().topLeftCorner(dimv_, dimv_) * dImD_dq_;
  MJtJinvImDCqv_().topLeftCorner(dimv_, dimv_).noalias() 
      += MJtJinv_().topRightCorner(dimv_, dimf_) * kkt_matrix.Cq();
  MJtJinvImDCqv_().topRightCorner(dimv_, dimv_).noalias() 
      = MJtJinv_().topRightCorner(dimv_, dimf_) * kkt_matrix.Cv();
  MJtJinvImDCqv_().bottomLeftCorner(dimf_, dimv_).noalias() 
      = MJtJinv_().bottomLeftCorner(dimf_, dimv_) * dImD_dq_;
  MJtJinvImDCqv_().bottomLeftCorner(dimf_, dimv_).noalias() 
      += MJtJinv_().bottomRightCorner(dimf_, dimf_) * kkt_matrix.Cq();
  MJtJinvImDCqv_().bottomRightCorner(dimf_, dimv_).noalias() 
      = MJtJinv_().bottomRightCorner(dimf_, dimf_) * kkt_matrix.Cv();
  MJtJinvImDC_().noalias() = MJtJinv_().leftCols(dimv_) * kkt_residual.dv_res;
  MJtJinvImDC_().noalias() += MJtJinv_().rightCols(dimf_) * kkt_residual.C();
  Qdvq_condensed_.noalias() = (- kkt_matrix.Qdvdv.diagonal()).asDiagonal() 
                                * MJtJinvImDCqv_().topLeftCorner(dimv_, dimv_);
  Qdvv_condensed_.noalias() = (- kkt_matrix.Qdvdv.diagonal()).asDiagonal() 
                                * MJtJinvImDCqv_().topRightCorner(dimv_, dimv_);
  Qfq_condensed_().noalias() = (- kkt_matrix.Qff().diagonal()).asDiagonal() 
                                * MJtJinvImDCqv_().bottomLeftCorner(dimf_, dimv_);
  Qfv_condensed_().noalias() = (- kkt_matrix.Qff().diagonal()).asDiagonal() 
                                * MJtJinvImDCqv_().bottomRightCorner(dimf_, dimv_);
  ldv_condensed_ = kkt_residual.ldv;
  ldv_condensed_.array() -= kkt_matrix.Qdvdv.diagonal().array() 
                              * MJtJinvImDC_().head(dimv_).array();
  lf_condensed_() = - kkt_residual.lf();
  lf_condensed_().array() -= kkt_matrix.Qff().diagonal().array() 
                              * MJtJinvImDC_().tail(dimf_).array();
  kkt_matrix.Qqq().noalias()
      -= MJtJinvImDCqv_().topLeftCorner(dimv_, dimv_).transpose() 
          * Qdvq_condensed_;
  kkt_matrix.Qqq().noalias()
      -= MJtJinvImDCqv_().bottomLeftCorner(dimf_, dimv_).transpose() 
          * Qfq_condensed_();
  kkt_matrix.Qqv().noalias()
      -= MJtJinvImDCqv_().topLeftCorner(dimv_, dimv_).transpose() 
          * Qdvv_condensed_;
  kkt_matrix.Qqv().noalias()
      -= MJtJinvImDCqv_().bottomLeftCorner(dimf_, dimv_).transpose() 
          * Qfv_condensed_();
  kkt_matrix.Qvq().noalias() = kkt_matrix.Qqv().transpose();
  kkt_matrix.Qvv().noalias()
      -= MJtJinvImDCqv_().topRightCorner(dimv_, dimv_).transpose() 
          * Qdvv_condensed_;
  kkt_matrix.Qvv().noalias()
      -= MJtJinvImDCqv_().bottomRightCorner(dimf_, dimv_).transpose() 
          * Qfv_condensed_();
  kkt_residual.lq().noalias()
      -= MJtJinvImDCqv_().topLeftCorner(dimv_, dimv_).transpose() 
          * ldv_condensed_;
  kkt_residual.lq().noalias()
      -= MJtJinvImDCqv_().bottomLeftCorner(dimf_, dimv_).transpose() 
          * lf_condensed_();
  kkt_residual.lv().noalias()
      -= MJtJinvImDCqv_().topRightCorner(dimv_, dimv_).transpose() 
          * ldv_condensed_;
  kkt_residual.lv().noalias()
      -= MJtJinvImDCqv_().bottomRightCorner(dimf_, dimv_).transpose() 
          * lf_condensed_();
  kkt_matrix.Fvq = - MJtJinvImDCqv_().topLeftCorner(dimv_, dimv_);
  kkt_matrix.Fvv = Eigen::MatrixXd::Identity(dimv_, dimv_) 
                    - MJtJinvImDCqv_().topRightCorner(dimv_, dimv_);
  kkt_residual.Fv().noalias() -= MJtJinvImDC_().head(dimv_);
}


inline void ImpulseDynamicsForwardEuler::computeCondensedDirection(
    const ImpulseKKTMatrix& kkt_matrix, const ImpulseKKTResidual& kkt_residual, 
    const SplitDirection& d_next, ImpulseSplitDirection& d) {
  d.ddv.noalias() = - MJtJinvImDCqv_().topRows(dimv_) * d.dx() 
                    - MJtJinvImDC_().head(dimv_);
  d.df().noalias() = - MJtJinvImDCqv_().bottomRows(dimf_) * d.dx()
                     - MJtJinvImDC_().tail(dimf_);
  ldv_condensed_.noalias() += d_next.dgmm();
  ldv_condensed_.noalias() += Qdvq_condensed_ * d.dq();
  ldv_condensed_.noalias() += Qdvv_condensed_ * d.dv();
  lf_condensed_().noalias() += Qfq_condensed_() * d.dq();
  lf_condensed_().noalias() += Qfv_condensed_() * d.dv();
  d.dbeta.noalias() = - MJtJinv_().topLeftCorner(dimv_, dimv_) * ldv_condensed_
                      - MJtJinv_().topRightCorner(dimv_, dimf_) * lf_condensed_();
  d.dmu().noalias() = - MJtJinv_().bottomLeftCorner(dimf_, dimv_) * ldv_condensed_
                      - MJtJinv_().bottomRightCorner(dimf_, dimf_) * lf_condensed_();
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
  robot.dRNEAPartialdFext(contact_status, dImD_df_full_);
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


inline Eigen::Block<Eigen::MatrixXd> ImpulseDynamicsForwardEuler::dImD_df_() {
  return dImD_df_full_.topLeftCorner(dimv_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseDynamicsForwardEuler::MJtJinv_() {
  return MJtJinv_full_.topLeftCorner(dimv_+dimf_, dimv_+dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsForwardEuler::MJtJinvImDCqv_() {
  return MJtJinvImDCqv_full_.topLeftCorner(dimv_+dimf_, 2*dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsForwardEuler::Qfq_condensed_() {
  return Qfq_condensed_full_.topLeftCorner(dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsForwardEuler::Qfv_condensed_() {
  return Qfv_condensed_full_.topLeftCorner(dimf_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsForwardEuler::MJtJinvImDC_() {
  return MJtJinvImDC_full_.head(dimv_+dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsForwardEuler::lf_condensed_() {
  return lf_condensed_full_.head(dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsForwardEuler::dImD_df_() const {
  return dImD_df_full_.topLeftCorner(dimv_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsForwardEuler::MJtJinv_() const {
  return MJtJinv_full_.topLeftCorner(dimv_+dimf_, dimv_+dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsForwardEuler::MJtJinvImDCqv_() const {
  return MJtJinvImDCqv_full_.topLeftCorner(dimv_+dimf_, 2*dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsForwardEuler::Qfq_condensed_() const {
  return Qfq_condensed_full_.topLeftCorner(dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsForwardEuler::Qfv_condensed_() const {
  return Qfv_condensed_full_.topLeftCorner(dimf_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseDynamicsForwardEuler::MJtJinvImDC_() const {
  return MJtJinvImDC_full_.head(dimv_+dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseDynamicsForwardEuler::lf_condensed_() const {
  return lf_condensed_full_.head(dimf_);
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HXX_ 