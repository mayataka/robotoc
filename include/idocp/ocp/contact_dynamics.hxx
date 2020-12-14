#ifndef IDOCP_CONTACT_DYNAMICS_HXX_
#define IDOCP_CONTACT_DYNAMICS_HXX_

#include "idocp/ocp/contact_dynamics.hpp"

#include <assert.h>

namespace idocp {

inline ContactDynamics::ContactDynamics(const Robot& robot) 
  : data_(robot),
    has_floating_base_(robot.has_floating_base()),
    has_active_contacts_(false) {
}


inline ContactDynamics::ContactDynamics() 
  : data_(),
    has_floating_base_(false),
    has_active_contacts_(false) {
}


inline ContactDynamics::~ContactDynamics() {
}


inline void ContactDynamics::linearizeContactDynamics(
    Robot& robot, const ContactStatus& contact_status, const double dtau, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) { 
  assert(dtau > 0);
  setContactStatus(contact_status);
  linearizeInverseDynamics(robot, contact_status, s, data_);
  linearizeContactConstraint(robot, contact_status, dtau, data_);
  // augment inverse dynamics constraint
  kkt_residual.lq().noalias() += dtau * data_.dIDdq().transpose() * s.beta;
  kkt_residual.lv().noalias() += dtau * data_.dIDdv().transpose() * s.beta;
  kkt_residual.la.noalias() += dtau * data_.dIDda.transpose() * s.beta;
  if (has_active_contacts_) {
    // We use an equivalence dIDdf_().transpose() = - dCda_(), to avoid
    // redundant calculation of dIDdf_().
    kkt_residual.lf().noalias() -= dtau * data_.dCda() * s.beta;
  }
  // augment floating base constraint
  if (has_floating_base_) {
    kkt_residual.lu().noalias() -= dtau * s.beta.tail(robot.dimu()); 
    kkt_residual.lu_passive.noalias() 
        -= dtau * s.beta.template head<kDimFloatingBase>(); 
    data_.u_passive = s.u_passive; 
    kkt_residual.lu_passive.noalias() += dtau * s.nu_passive;
  }
  else {
    kkt_residual.lu().noalias() -= dtau * s.beta; 
  }
  // augment contact constraint
  if (has_active_contacts_) {
    kkt_residual.lq().noalias() 
        += dtau * data_.dCdq().transpose() * s.mu_stack();
    kkt_residual.lv().noalias() 
        += dtau * data_.dCdv().transpose() * s.mu_stack();
    kkt_residual.la.noalias() += dtau * data_.dCda().transpose() * s.mu_stack();
  }
}


inline void ContactDynamics::linearizeInverseDynamics(
    Robot& robot, const ContactStatus& contact_status, 
    const SplitSolution& s, ContactDynamicsData& data) {
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, data.ID_full());
  if (robot.has_floating_base()) {
    data.ID_passive().noalias() -= s.u_passive;
    data.ID().noalias() -= s.u;
  }
  else {
    data.ID_full().noalias() -= s.u;
  }
  robot.RNEADerivatives(s.q, s.v, s.a, data.dIDdq(), data.dIDdv(), data.dIDda);
}


inline void ContactDynamics::linearizeContactConstraint(
    Robot& robot, const ContactStatus& contact_status, const double dtau, 
    ContactDynamicsData& data) {
  assert(dtau > 0);
  robot.computeBaumgarteResidual(contact_status, dtau, 
                                 contact_status.contactPoints(), data.C());
  robot.computeBaumgarteDerivatives(contact_status, dtau, data.dCdq(), 
                                    data.dCdv(), data.dCda());
}


inline void ContactDynamics::condenseContactDynamics(
    Robot& robot, const ContactStatus& contact_status, const double dtau,
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) {
  assert(dtau > 0);
  robot.computeMJtJinv(data_.dIDda, data_.dCda(), data_.MJtJinv());
  condensing(robot, dtau, data_, kkt_matrix, kkt_residual);
}


inline void ContactDynamics::condensing(const Robot& robot, const double dtau, 
                                        ContactDynamicsData& data, 
                                        SplitKKTMatrix& kkt_matrix, 
                                        SplitKKTResidual& kkt_residual) {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dim_passive = robot.dim_passive();
  data.MJtJinv_dIDCdqv().noalias() = data.MJtJinv() * data.dIDCdqv();
  data.MJtJinv_IDC().noalias() = data.MJtJinv() * data.IDC();
  data.Qafqv().noalias() 
      = (- kkt_matrix.Qaaff().diagonal()).asDiagonal() * data.MJtJinv_dIDCdqv();
  data.Qafu_full().noalias() 
      = kkt_matrix.Qaaff().diagonal().asDiagonal() 
          * data.MJtJinv().leftCols(dimv);
  data.la() = kkt_residual.la;
  data.lf() = - kkt_residual.lf();
  data.laf().noalias() 
      -= kkt_matrix.Qaaff().diagonal().asDiagonal() * data.MJtJinv_IDC();
  kkt_matrix.Qxx().noalias() 
      -= data.MJtJinv_dIDCdqv().transpose() * data.Qafqv();
  kkt_matrix.Qxu_full().noalias() 
      -= data.MJtJinv_dIDCdqv().transpose() * data.Qafu_full();
  kkt_residual.lx().noalias() 
      -= data.MJtJinv_dIDCdqv().transpose() * data.laf();
  if (robot.has_floating_base()) {
    kkt_matrix.Quu_full().noalias() 
        += data.MJtJinv().topRows(dimv) * data.Qafu_full();
    kkt_residual.lu_passive.noalias() 
        += data.MJtJinv().template topRows<kDimFloatingBase>() * data.laf();
    kkt_residual.lu().noalias() 
        += data.MJtJinv().middleRows(kDimFloatingBase, dimu) * data.laf();
    kkt_residual.lu().noalias() 
        -= kkt_matrix.Quu_passive_bottomLeft() * data.u_passive;
    kkt_residual.lx().noalias() -= kkt_matrix.Qxu_passive() * data.u_passive;
  }
  else {
    kkt_matrix.Quu().noalias() += data.MJtJinv().topRows(dimv) * data.Qafu();
    kkt_residual.lu().noalias() += data.MJtJinv().topRows(dimv) * data.laf();
  }
  kkt_matrix.Fvq() = - dtau * data.MJtJinv_dIDCdqv().topLeftCorner(dimv, dimv);
  kkt_matrix.Fvv().noalias() 
      = Eigen::MatrixXd::Identity(dimv, dimv) 
          - dtau * data.MJtJinv_dIDCdqv().topRightCorner(dimv, dimv);
  if (robot.has_floating_base()) {
    kkt_matrix.Fvu() = dtau * data.MJtJinv().block(0, dim_passive, dimv, dimu);
    kkt_residual.Fv().noalias() 
        -= dtau * data.MJtJinv().topLeftCorner(dimv, dim_passive) 
                * data.u_passive;
    kkt_residual.Fv().noalias() -= dtau * data.MJtJinv_IDC().head(dimv);
  }
  else {
    kkt_matrix.Fvu() = dtau * data.MJtJinv().topLeftCorner(dimv, dimv);
    kkt_residual.Fv().noalias() -= dtau * data.MJtJinv_IDC().head(dimv);
  }
}


inline void ContactDynamics::computeCondensedPrimalDirection(
    const Robot& robot, SplitDirection& d) {
  expansionPrimal(robot, data_, d);
}


template <typename VectorType>
inline void ContactDynamics::computeCondensedDualDirection(
    const Robot& robot, const double dtau, const SplitKKTMatrix& kkt_matrix, 
    const SplitKKTResidual& kkt_residual, 
    const Eigen::MatrixBase<VectorType>& dgmm, SplitDirection& d) {
  assert(dtau > 0);
  assert(dgmm.size() == robot.dimv());
  expansionDual(robot, dtau, data_, kkt_matrix, kkt_residual, dgmm, d);
}


inline void ContactDynamics::expansionPrimal(const Robot& robot, 
                                             const ContactDynamicsData& data, 
                                             SplitDirection& d) {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  if (robot.has_floating_base()) {
    d.du_passive = - data.u_passive;
    d.daf().noalias() = - data.MJtJinv_dIDCdqv() * d.dx();
    d.daf().noalias() 
        += data.MJtJinv().template leftCols<kDimFloatingBase>() * d.du_passive;
    d.daf().noalias() 
        += data.MJtJinv().middleCols(kDimFloatingBase, dimu) * d.du();
    d.daf().noalias() -= data.MJtJinv_IDC();
  }
  else {
    d.daf().noalias() = - data.MJtJinv_dIDCdqv() * d.dx();
    d.daf().noalias() += data.MJtJinv().leftCols(dimv) * d.du();
    d.daf().noalias() -= data.MJtJinv_IDC();
  }
  d.df().array() *= -1;
}


template <typename VectorType>
inline void ContactDynamics::expansionDual(
    const Robot& robot, const double dtau, ContactDynamicsData& data, 
    const SplitKKTMatrix& kkt_matrix, const SplitKKTResidual& kkt_residual, 
    const Eigen::MatrixBase<VectorType>& dgmm, SplitDirection& d) {
  assert(dtau > 0);
  assert(dgmm.size() == robot.dimv());
  const int dimv = robot.dimv();
  if (robot.has_floating_base()) {
    d.dnu_passive = kkt_residual.lu_passive;
    d.dnu_passive.noalias() += kkt_matrix.Quu_passive_topLeft() * d.du_passive;
    d.dnu_passive.noalias() += kkt_matrix.Quu_passive_topRight() * d.du();
    d.dnu_passive.noalias() += kkt_matrix.Qxu_passive().transpose() * d.dx();
    d.dnu_passive.noalias() 
        += dtau * data.MJtJinv().topLeftCorner(kDimFloatingBase, dimv) * dgmm;
    d.dnu_passive.array() *= - (1/dtau);
  }
  data.laf().noalias() += data.Qafqv() * d.dx();
  if (robot.has_floating_base()) {
    data.laf().noalias() += data.Qafu_passive() * d.du_passive;
    data.laf().noalias() += data.Qafu() * d.du();
  }
  else {
    data.laf().noalias() += data.Qafu() * d.du();
  }
  data.la().noalias() += dtau * dgmm;
  d.dbetamu().noalias() = - data.MJtJinv() * data.laf() * (1/dtau);
}


inline void ContactDynamics::computeContactDynamicsResidual(
    Robot& robot, const ContactStatus& contact_status, const double dtau,
    const SplitSolution& s) {
  setContactStatus(contact_status);
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, data_.ID_full());
  if (has_floating_base_) {
    data_.ID_passive().noalias() -= s.u_passive;
    data_.ID().noalias() -= s.u;
    data_.u_passive = s.u_passive;
  }
  else {
    data_.ID_full().noalias() -= s.u;
  }
  robot.computeBaumgarteResidual(contact_status, dtau, 
                                 contact_status.contactPoints(), data_.C());
}


inline double ContactDynamics::l1NormContactDynamicsResidual(
    const double dtau) const {
  assert(dtau > 0);
  if (has_floating_base_) {
    return (dtau * (data_.IDC().lpNorm<1>() + data_.u_passive.lpNorm<1>()));
  }
  else {
    return (dtau * data_.IDC().lpNorm<1>());
  }
}


inline double ContactDynamics::squaredNormContactDynamicsResidual(
    const double dtau) const {
  assert(dtau > 0);
  if (has_floating_base_) {
    return (dtau * dtau * (data_.IDC().squaredNorm() 
                            + data_.u_passive.squaredNorm()));
  }
  else {
    return (dtau * dtau * data_.IDC().squaredNorm());
  }
}


inline void ContactDynamics::setContactStatus(
    const ContactStatus& contact_status) {
  data_.setContactStatus(contact_status);
  has_active_contacts_ = contact_status.hasActiveContacts();
}

} // namespace idocp 

#endif // IDOCP_CONTACT_DYNAMICS_HXX_ 