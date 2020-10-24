#ifndef IDOCP_CONTACT_DYNAMICS_HXX_
#define IDOCP_CONTACT_DYNAMICS_HXX_

#include "idocp/ocp/contact_dynamics.hpp"

#include <assert.h>

namespace idocp {

inline ContactDynamicsForwardEuler::ContactDynamicsForwardEuler(
    const Robot& robot) 
  : dIDCdqv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                        2*robot.dimv())),
    MJtJinv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                        robot.dimv()+robot.max_dimf())), 
    MJtJinvImDCqv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                              2*robot.dimv())), 
    Qafqv_condensed_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                                robot.dimv()+robot.max_dimf())), 
    Qafu_condensed_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                               robot.dimv())), 
    MJtJinv_IDC_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    laf_condensed_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    dimv_(robot.dimv()),
    dimf_(0) {
}


inline ContactDynamicsForwardEuler::ContactDynamicsForwardEuler() 
  : dIDCdqv_full_(),
    MJtJinv_full_(), 
    MJtJinvImDCqv_full_(), 
    Qafqv_condensed_full_(), 
    Qafu_condensed_full_(), 
    MJtJinv_IDC_full_(),
    laf_condensed_full_(),
    dimv_(0),
    dimf_(0) {
}


inline ContactDynamicsForwardEuler::~ContactDynamicsForwardEuler() {
}


inline void ContactDynamicsForwardEuler::linearizeContactDynamics(
    Robot& robot, const ContactStatus& contact_status, const double dtau, 
    const ImpulseSplitSolution& s, ImpulseKKTMatrix& kkt_matrix, 
    ImpulseKKTResidual& kkt_residual) { 
  setContactStatus(contact_status);
  linearizeInverseDynamics(robot, contact_status, s, kkt_residual);
  linearizeContactConstraint(robot, contact_status, dtau, kkt_matrix, 
                              kkt_residual);
  // augment inverse dynamics constraint
  kkt_residual.la().noalias() += dtau * dIDda_().transpose() * s.beta;
  if (has_active_contacts_) {
    kkt_residual.lf().noalias() -= dtau * dCda_() * s.beta;
  }
  kkt_residual.lq().noalias() += dtau * dIDdq_().transpose() * s.beta;
  kkt_residual.lv().noalias() += dtau * dIDdv_().transpose() * s.beta;
  kkt_residual.lu.noalias() -= dtau * s.beta; 
  // augment floating base constraint
  if (has_floating_base_) {
    computeFloatingBaseConstraintResidual(robot, dtau, s, kkt_residual);
    kkt_residual.lu.template head<kDimFloatingBase>().noalias() 
        += dtau * s.mu_floating_base();
  }
  // augment contact constraint
  if (has_active_contacts_) {
    kkt_residual.la().noalias() += dCda_().transpose() * s.mu_contacts();
    kkt_residual.lq().noalias() += dCdq_().transpose() * s.mu_contacts();
    kkt_residual.lv().noalias() += dCdv_().transpose() * s.mu_contacts();
  }
}


inline void ContactDynamicsForwardEuler::condenseContactDynamics(
    Robot& robot, const ContactStatus& contact_status, const double dtau,
    const ImpulseSplitSolution& s, ImpulseKKTMatrix& kkt_matrix, 
    ImpulseKKTResidual& kkt_residual) {
  setContactStatus(contact_status);
  linearizeInverseDynamics(robot, contact_status, s, kkt_residual);
  linearizeContactConstraint(robot, contact_status, dtau, kkt_matrix, 
                              kkt_residual);
  robot.computeMJtJinv(contact_status, dIDda_(), dCda(), MJtJinv_());
  MJtJinv_dIDCdqv_().noalias() = MJtJinv_() * dIDCdqv_();
  MJtJinv_IDC_().noalias() = MJtJinv_().leftCols(dimv_) * kkt_residual.u_res;
  MJtJinv_IDC_().noalias() += MJtJinv_().rightCols(dimf_) * kkt_residual.C();
  Qaq_condensed_().noalias() = (- kkt_matrix.Qaa().diagonal()).asDiagonal() 
                                * MJtJinv_dIDCdqv_().topLeftCorner(dimv_, dimv_);
  Qav_condensed_().noalias() = (- kkt_matrix.Qaa().diagonal()).asDiagonal() 
                                * MJtJinv_dIDCdqv_().topRightCorner(dimv_, dimv_);
  Qfq_condensed_().noalias() = (- kkt_matrix.Qff().diagonal()).asDiagonal() 
                                * MJtJinv_dIDCdqv_().bottomLeftCorner(dimf_, dimv_);
  Qfv_condensed_().noalias() = (- kkt_matrix.Qff().diagonal()).asDiagonal() 
                                * MJtJinv_dIDCdqv_().bottomRightCorner(dimf_, dimv_);
  Qau_condensed_().noalias() = (- kkt_matrix.Qaa().diagonal()).asDiagonal() 
                                * MJtJinv_().topLeftCorner(dimv_, dimv_);
  Qfu_condensed_().noalias() = (- kkt_matrix.Qff().diagonal()).asDiagonal() 
                                * MJtJinv_().bottomLeftCorner(dimv_, dimv_);
  la_condensed_() = kkt_residual.la();
  la_condensed_().array() -= kkt_matrix.Qaa().diagonal().array() 
                              * MJtJinv_IDC_().head(dimv_).array();
  lf_condensed_() = - kkt_residual.lf();
  lf_condensed_().array() -= kkt_matrix.Qff().diagonal().array() 
                              * MJtJinv_IDC_().tail(dimf_).array();
  kkt_matrix.Qxx().noalias() -= MJtJinv_dIDCdqv_().transpose() * Qafqv_condensed_();
  kkt_matrix.Qxu().noalias() -= MJtJinv_dIDCdqv_().transpose() * Qafu_condensed_();
  kkt_residual.lx().noalias() -= MJtJinv_dIDCdqv_().transpose() * laf_condensed_();
  kkt_matrix.Quu().noalias() += MJtJinv_().topRows(dimv_).transpose() * Qafu_condensed_();
  kkt_residual.lu().noalias() += MJtJinv_().topRows(dimv_).transpose() * laf_condensed_();
  if (has_floating_base_) {
    kkt_residual.lu_L().noalias() 
        -= kkt_matrix.Quu_LU() * s.u.template head<kDimFloatingBase>();
  }
  kkt_matrix.Fvq = - MJtJinv_dIDCdqv_().topLeftCorner(dimv_, dimv_);
  kkt_matrix.Fvv = dtau * Eigen::MatrixXd::Identity(dimv_, dimv_) 
                    - MJtJinv_dIDCdqv_().topRightCorner(dimv_, dimv_);
  kkt_residual.Fv().noalias() -= dtau * MJtJinv_IDC_().head(dimv_);
  kkt_residual.Fv().noalias() 
      -= dtau * MJtJinv_().topLeftCorner(dimv_, dim_passive_) 
              * s.u.head(dim_passive_);
}


inline void ContactDynamicsForwardEuler::computeCondensedDirection(
    const double dtau, const ImpulseKKTMatrix& kkt_matrix, 
    const ImpulseKKTResidual& kkt_residual, const SplitDirection& d_next, 
    ImpulseSplitDirection& d) {
  d.du.template head<kDimFloatingBase>() 
      = - s.u.template head<kDimFloatingBase>();
  d.dmu_floating_base().noalias()
      = (kkt_residual.lu_U() 
          + kkt_matrix.Quu_UU() * d.du.template head<kDimFloatingBase>()
          + kkt_matrix.Quu_UL() * d.du.tail(dimv_-dim_passive_)
          + kkt_matrix.Qua_U() * d.dq() + kkt_matrix.Quv_U() * d.dv() 
          + dtau * MJtJinv_().topLeftCorner(dim_passive_, dimv_) * d_next.dgmm()
          ) / dtau;
  d.da().noalias() = - MJtJinv_dIDCdqv_().topRows(dimv_) * d.dx() 
                     + MJtJinv_().topLeftCorner(dimv_, dimv_) * d.du
                     - MJtJinv_IDC_().head(dimv_);
  d.df().noalias() = MJtJinv_dIDCdqv_().bottomRows(dimf_) * d.dx()
                     - MJtJinv_().bottomLeftCorner(dimf_, dimv_) * d.du
                     + MJtJinv_IDC_().tail(dimf_);
  laf_condensed_().noalias() += Qafqv_condensed_() * d.dx();
  laf_condensed_().noalias() += Qafu_condensed_() * d.du;
  la_condensed_().noalias() += dtau * d_next.dgmm();
  d.dbeta.noalias() = - MJtJinv_().topRows(dimv_) * ldvf_condensed_();
  d.dmu().noalias() = - MJtJinv_().bottomRows(dimf_) * ldvf_condensed_();
}


inline void ContactDynamicsForwardEuler::computeContactDynamicsResidual(
    Robot& robot, const ContactStatus& contact_status, const double dtau,
    const ImpulseSplitSolution& s, ImpulseKKTResidual& kkt_residual) {
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual.u_res);
  kkt_residual.u_res.noalias() -= s.u;
  robot.computeContactyResidual(contact_status, kkt_residual.C());
}


inline void ContactDynamicsForwardEuler::linearizeInverseContactDynamics(
    Robot& robot, const ContactStatus& contact_status, 
    const ImpulseSplitSolution& s, ImpulseKKTResidual& kkt_residual) {
  setContactForces(robot, contact_status, s);
  computeInverseDynamicsResidual(robot, s, kkt_residual);
  robot.RNEADerivatives(s.q, s.v, s.a, dID_dq_(), dID_dv_(), dID_da_());
}


inline void ContactDynamicsForwardEuler::linearizeContactConstraint(
    Robot& robot, const ContactStatus& contact_status, 
    ImpulseKKTMatrix& kkt_matrix, ImpulseKKTResidual& kkt_residual) {
  assert(dtau > 0);
  computeContactConstraintResidual(robot, contact_status, dtau, kkt_residual);
  robot.computeBaumgarteDerivatives(contact_status, dtau, dtau, 
                                    Cq_(), Cv_(), Ca_());
}


inline double ContactDynamicsForwardEuler::l1NormContactDynamicsResidual(
    const ImpulseKKTResidual& kkt_residual) {
  return (kkt_residual.u_res.lpNorm<1>() + kkt_residual.C().lpNorm<1>());
}


inline double ContactDynamicsForwardEuler::squaredNormContactDynamicsResidual(
    const ImpulseKKTResidual& kkt_residual) {
  return (kkt_residual.dv_res.squaredNorm() + kkt_residual.C().squaredNorm());
}


inline void ContactDynamicsForwardEuler::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamicsForwardEuler::dIDCdqv_() {
  return dIDCdqv_full_.topLeftCorner(dimv_+dimf_, 2*dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamicsForwardEuler::dIDdq_() {
  return dIDCdqv_full_.topLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamicsForwardEuler::dIDdv_() {
  return dIDCdqv_full_.block(0, dimv_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamicsForwardEuler::dCdq_() {
  return dIDCdqv_full_.block(dimv_, 0, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamicsForwardEuler::dCdv_() {
  return dIDCdqv_full_.block(dimv_, dimv_, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamicsForwardEuler::dCda_() {
  return dIDCdqv_full_.topLeftCorner(dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamicsForwardEuler::MJtJinv_() {
  return MJtJinv_full_.topLeftCorner(dimv_+dimf_, dimv_+dimf_);
}

inline Eigen::Block<Eigen::MatrixXd> 
ContactDynamicsForwardEuler::MJtJinv_dIDCdqv_() {
  return MJtJinv_dIDCdqv_full_.topLeftCorner(dimv_+dimf_, 2*dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ContactDynamicsForwardEuler::Qafqv_condensed_() {
  return Qafqv_condensed_full_.topLeftCorner(dimv_+dimf_, 2*dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ContactDynamicsForwardEuler::Qafu_condensed_() {
  return Qafqv_condensed_full_.topLeftCorner(dimv_+dimf_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ContactDynamicsForwardEuler::MJtJinv_IDC_() {
  return MJtJinv_IDC_full_.head(dimv_+dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ContactDynamicsForwardEuler::laf_condensed_() {
  return laf_condensed_full_.head(dimv_+dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ContactDynamicsForwardEuler::la_condensed_() {
  return laf_condensed_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ContactDynamicsForwardEuler::lf_condensed_() {
  return laf_condensed_full_.segment(dimv_, dimf_);
}

} // namespace idocp 

#endif // IDOCP_CONTACT_DYNAMICS_HXX_ 