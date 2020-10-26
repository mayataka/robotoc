#ifndef IDOCP_CONTACT_DYNAMICS_HXX_
#define IDOCP_CONTACT_DYNAMICS_HXX_

#include "idocp/ocp/contact_dynamics.hpp"

#include <assert.h>

namespace idocp {

inline ContactDynamics::ContactDynamics(const Robot& robot) 
  : dIDda_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dCda_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
    dIDCdqv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                        2*robot.dimv())),
    MJtJinv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                        robot.dimv()+robot.max_dimf())), 
    MJtJinv_dIDCdqv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                                2*robot.dimv())), 
    Qafqv_condensed_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                                2*robot.dimv())), 
    Qafu_condensed_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                               robot.dimv())), 
    u_passive_(Eigen::VectorXd::Zero(robot.dim_passive())),
    IDC_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    MJtJinv_IDC_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    laf_condensed_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    dimv_(robot.dimv()),
    dimu_(robot.dimu()),
    dimf_(0),
    dim_passive_(robot.dim_passive()),
    has_floating_base_(robot.has_floating_base()),
    has_active_contacts_(false) {
}


inline ContactDynamics::ContactDynamics() 
  : dIDda_(),
    dCda_full_(),
    dIDCdqv_full_(),
    MJtJinv_full_(), 
    MJtJinv_dIDCdqv_full_(), 
    Qafqv_condensed_full_(), 
    Qafu_condensed_full_(), 
    u_passive_(),
    IDC_full_(),
    MJtJinv_IDC_full_(),
    laf_condensed_full_(),
    dimv_(0),
    dimu_(0),
    dimf_(0),
    dim_passive_(0),
    has_floating_base_(false),
    has_active_contacts_(false) {
}


inline ContactDynamics::~ContactDynamics() {
}


inline void ContactDynamics::linearizeContactDynamics(
    Robot& robot, const ContactStatus& contact_status, const double dtau, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) { 
  assert(dtau > 0);
  setContactStatus(contact_status);
  linearizeInverseDynamics(robot, contact_status, s, kkt_residual);
  linearizeContactConstraint(robot, contact_status, dtau, kkt_matrix, 
                              kkt_residual);
  // augment inverse dynamics constraint
  kkt_residual.lq().noalias() += dtau * dIDdq_().transpose() * s.beta;
  kkt_residual.lv().noalias() += dtau * dIDdv_().transpose() * s.beta;
  kkt_residual.la.noalias() += dtau * dIDda_.transpose() * s.beta;
  if (has_active_contacts_) {
    // We use an equivalence, dIDdf_().transpose() = - dCda_(), to avoid
    // redundant calculation of dIDdf_().
    kkt_residual.lf().noalias() -= dtau * dCda_() * s.beta;
  }
  // augment floating base constraint
  if (has_floating_base_) {
    kkt_residual.lu().noalias() -= dtau * s.beta.tail(dimu_); 
    kkt_residual.lu_passive.noalias() -= dtau * s.beta.template head<kDimFloatingBase>(); 
    u_passive_ = s.u_passive; 
    kkt_residual.lu_passive.noalias() += dtau * s.nu_passive;
  }
  else {
    kkt_residual.lu().noalias() -= dtau * s.beta; 
  }
  // augment contact constraint
  if (has_active_contacts_) {
    kkt_residual.lq().noalias() += dtau * dCdq_().transpose() * s.mu_stack();
    kkt_residual.lv().noalias() += dtau * dCdv_().transpose() * s.mu_stack();
    kkt_residual.la.noalias() += dtau * dCda_().transpose() * s.mu_stack();
  }
}


inline void ContactDynamics::condenseContactDynamics(
    Robot& robot, const ContactStatus& contact_status, const double dtau,
    const SplitSolution& s, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) {
  assert(dtau > 0);
  linearizeContactDynamics(robot, contact_status, dtau, s, 
                           kkt_matrix, kkt_residual);
  robot.computeMJtJinv(contact_status, dIDda_, dCda_(), MJtJinv_());
  MJtJinv_dIDCdqv_().noalias() = MJtJinv_() * dIDCdqv_();
  MJtJinv_IDC_().noalias() = MJtJinv_() * IDC_();
  Qafqv_condensed_().noalias() = (- kkt_matrix.Qaaff().diagonal()).asDiagonal() * MJtJinv_dIDCdqv_();
  Qafu_condensed_().noalias() = (- kkt_matrix.Qaaff().diagonal()).asDiagonal() * MJtJinv_().leftCols(dimv_);
  la_condensed_() = kkt_residual.la;
  lf_condensed_() = - kkt_residual.lf();
  la_condensed_().noalias() -= kkt_matrix.Qaaff().diagonal().asDiagonal() * MJtJinv_IDC_();
  kkt_matrix.Qxx().noalias() -= MJtJinv_dIDCdqv_().transpose() * Qafqv_condensed_();
  kkt_matrix.Qxu_full().noalias() -= MJtJinv_dIDCdqv_().transpose() * Qafu_condensed_();
  kkt_residual.lx().noalias() -= MJtJinv_dIDCdqv_().transpose() * laf_condensed_();
  if (has_floating_base_) {
    kkt_matrix.Quu_full() = MJtJinv_().topRows(dimv_).transpose() * Qafu_condensed_();
    kkt_residual.lu_passive.noalias() 
        += MJtJinv_().topRows(dimv_).transpose().template topRows<kDimFloatingBase>() * laf_condensed_();
    kkt_residual.lu_passive.noalias() 
        += MJtJinv_().topRows(dimv_).transpose().bottomRows(dimu_) * laf_condensed_();
    kkt_residual.lu().noalias() -= kkt_matrix.Quu_passive_drive().transpose() * s.u_passive;
  }
  else {
    kkt_matrix.Quu().noalias() += MJtJinv_().topRows(dimv_).transpose() * Qafu_condensed_();
    kkt_residual.lu().noalias() += MJtJinv_().topRows(dimv_).transpose() * laf_condensed_();
  }
  kkt_matrix.Fvq() = - dtau * MJtJinv_dIDCdqv_().topLeftCorner(dimv_, dimv_);
  kkt_matrix.Fvv() = Eigen::MatrixXd::Identity(dimv_, dimv_) 
                      - dtau * MJtJinv_dIDCdqv_().topRightCorner(dimv_, dimv_);
  kkt_matrix.Fvu() = dtau * MJtJinv_().topLeftCorner(dimv_, dimv_).topRightCorner(dimv_, dimu_);
  kkt_residual.Fv().noalias() -= dtau * MJtJinv_IDC_().head(dimv_);
  kkt_residual.Fv().noalias() -= dtau * MJtJinv_().topLeftCorner(dimv_, dim_passive_) * s.u_passive;
}


template <typename VectorType>
inline void ContactDynamics::computeCondensedDirection(
    const double dtau, const KKTMatrix& kkt_matrix, 
    const KKTResidual& kkt_residual, const Eigen::MatrixBase<VectorType>& dgmm, 
    SplitDirection& d) {
  assert(dgmm.size() == dimv_);
  if (has_floating_base_) {
    d.du_passive = - u_passive_;
    d.dnu_passive.noalias()
        = kkt_residual.lu_passive + kkt_matrix.Quu_passive() * d.du_passive
            + kkt_matrix.Quu_passive_drive() * d.du()
            + kkt_matrix.Quq_passive() * d.dq() 
            + kkt_matrix.Quv_passive() * d.dv() 
            + dtau * MJtJinv_().topLeftCorner(dim_passive_, dimv_) * dgmm;
    d.dnu_passive.array() /= - dtau;
    d.daf().noalias() = - MJtJinv_dIDCdqv_() * d.dx() 
                        + MJtJinv_().topLeftCorner(6, dimv_) * d.du_passive
                        + MJtJinv_().bottomLeftCorner(dimu_, dimv_) * d.du()
                        - MJtJinv_IDC_();
  }
  else {
    d.daf().noalias() = - MJtJinv_dIDCdqv_() * d.dx() 
                        + MJtJinv_().leftCols(dimv_) * d.du() - MJtJinv_IDC_();
  }
  laf_condensed_().noalias() += Qafqv_condensed_() * d.dx();
  laf_condensed_().noalias() += Qafu_condensed_() * d.du();
  la_condensed_().noalias() += dtau * dgmm;
  d.dbetamu().noalias() = - MJtJinv_() * laf_condensed_() / dtau;
}


inline void ContactDynamics::computeContactDynamicsResidual(
    Robot& robot, const ContactStatus& contact_status, const double dtau,
    const SplitSolution& s, KKTResidual& kkt_residual) {
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, ID_());
  if (has_floating_base_) {
    ID_().template head<kDimFloatingBase>().noalias() -= s.u_passive;
    ID_().tail(dimu_).noalias() -= s.u;
  }
  else {
    ID_().noalias() -= s.u;
  }
  robot.computeBaumgarteResidual(contact_status, dtau, C_());
}


inline void ContactDynamics::linearizeInverseDynamics(
    Robot& robot, const ContactStatus& contact_status, 
    const SplitSolution& s, KKTResidual& kkt_residual) {
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, ID_());
  if (has_floating_base_) {
    ID_().template head<kDimFloatingBase>().noalias() -= s.u_passive;
    ID_().tail(dimu_).noalias() -= s.u;
  }
  else {
    ID_().noalias() -= s.u;
  }
  robot.RNEADerivatives(s.q, s.v, s.a, dIDdq_(), dIDdv_(), dIDda_);
}


inline void ContactDynamics::linearizeContactConstraint(
    Robot& robot, const ContactStatus& contact_status, const double dtau,
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  assert(dtau > 0);
  robot.computeBaumgarteResidual(contact_status, dtau, C_());
  robot.computeBaumgarteDerivatives(contact_status, dtau, dCdq_(), dCdv_(), dCda_());
}


inline double ContactDynamics::l1NormContactDynamicsResidual(
    const double dtau) const {
  assert(dtau > 0);
  return (dtau * (IDC_().lpNorm<1>() + u_passive_.lpNorm<1>()));
}


inline double ContactDynamics::squaredNormContactDynamicsResidual(
    const double dtau) const {
  assert(dtau > 0);
  return (dtau * dtau * (IDC_().squaredNorm() + u_passive_.squaredNorm()));
}


inline void ContactDynamics::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
  has_active_contacts_ = contact_status.hasActiveContacts();
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamics::dIDCdqv_() {
  return dIDCdqv_full_.topLeftCorner(dimv_+dimf_, 2*dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamics::dIDdq_() {
  return dIDCdqv_full_.topLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamics::dIDdv_() {
  return dIDCdqv_full_.block(0, dimv_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamics::dCdq_() {
  return dIDCdqv_full_.block(dimv_, 0, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamics::dCdv_() {
  return dIDCdqv_full_.block(dimv_, dimv_, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamics::dCda_() {
  return dCda_full_.topLeftCorner(dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamics::MJtJinv_() {
  return MJtJinv_full_.topLeftCorner(dimv_+dimf_, dimv_+dimf_);
}

inline Eigen::Block<Eigen::MatrixXd> ContactDynamics::MJtJinv_dIDCdqv_() {
  return MJtJinv_dIDCdqv_full_.topLeftCorner(dimv_+dimf_, 2*dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamics::Qafqv_condensed_() {
  return Qafqv_condensed_full_.topLeftCorner(dimv_+dimf_, 2*dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ContactDynamics::Qafu_condensed_() {
  return Qafqv_condensed_full_.topLeftCorner(dimv_+dimf_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ContactDynamics::IDC_() {
  return IDC_full_.head(dimv_+dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> ContactDynamics::IDC_() const {
  return IDC_full_.head(dimv_+dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ContactDynamics::ID_() {
  return IDC_full_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ContactDynamics::ID_() const {
  return IDC_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ContactDynamics::C_() {
  return IDC_full_.segment(dimv_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ContactDynamics::C_() const {
  return IDC_full_.segment(dimv_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ContactDynamics::MJtJinv_IDC_() {
  return MJtJinv_IDC_full_.head(dimv_+dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ContactDynamics::laf_condensed_() {
  return laf_condensed_full_.head(dimv_+dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ContactDynamics::la_condensed_() {
  return laf_condensed_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ContactDynamics::lf_condensed_() {
  return laf_condensed_full_.segment(dimv_, dimf_);
}

} // namespace idocp 

#endif // IDOCP_CONTACT_DYNAMICS_HXX_ 