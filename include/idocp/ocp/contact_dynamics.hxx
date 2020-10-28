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
  kkt_residual.lq().noalias() += dtau * data_.dIDdq().transpose() * s.beta;
  kkt_residual.lv().noalias() += dtau * data_.dIDdv().transpose() * s.beta;
  kkt_residual.la.noalias() += dtau * data_.dIDda.transpose() * s.beta;
  if (has_active_contacts_) {
    // We use an equivalence, dIDdf_().transpose() = - dCda_(), to avoid
    // redundant calculation of dIDdf_().
    kkt_residual.lf().noalias() -= dtau * data_.dCda() * s.beta;
  }
  // augment floating base constraint
  if (has_floating_base_) {
    kkt_residual.lu().noalias() -= dtau * s.beta.tail(dimu_); 
    kkt_residual.lu_passive.noalias() -= dtau * s.beta.template head<kDimFloatingBase>(); 
    data_.u_passive_ = s.u_passive; 
    kkt_residual.lu_passive.noalias() += dtau * s.nu_passive;
  }
  else {
    kkt_residual.lu().noalias() -= dtau * s.beta; 
  }
  // augment contact constraint
  if (has_active_contacts_) {
    kkt_residual.lq().noalias() += dtau * data_.dCdq().transpose() * s.mu_stack();
    kkt_residual.lv().noalias() += dtau * data_.dCdv().transpose() * s.mu_stack();
    kkt_residual.la.noalias() += dtau * data_.dCda().transpose() * s.mu_stack();
  }
}


inline void ContactDynamics::condenseContactDynamics(
    Robot& robot, const ContactStatus& contact_status, const double dtau,
    const SplitSolution& s, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) {
  assert(dtau > 0);
  linearizeContactDynamics(robot, contact_status, dtau, s, 
                           kkt_matrix, kkt_residual);
  robot.computeMJtJinv(contact_status, data_.dIDda, data_.dCda(), data_.MJtJinv());
  condensing(dtau, data_, kkt_matrix, kkt_residual);
}


inline void ContactDynamics::condensing(const double dtau, 
                                        ContactDynamicsData& data,
                                        KKTMatrix& kkt_matrix, 
                                        KKTResidual& kkt_residual) const {
  data.MJtJinv_dIDCdqv().noalias() = data.MJtJinv() * data.dIDCdqv();
  data.MJtJinv_IDC().noalias() = data.MJtJinv() * data.IDC();
  data.Qafqv().noalias() = (- kkt_matrix.Qaaff().diagonal()).asDiagonal() * data.MJtJinv_dIDCdqv();
  data.Qafu_full().noalias() = (- kkt_matrix.Qaaff().diagonal()).asDiagonal() * data.MJtJinv().leftCols(dimv_);
  data.la() = kkt_residual.la;
  data.lf() = - kkt_residual.lf();
  data.laf().noalias() -= kkt_matrix.Qaaff().diagonal().asDiagonal() * data.MJtJinv_IDC();
  kkt_matrix.Qxx().noalias() -= data.MJtJinv_dIDCdqv().transpose() * data.Qafqv();
  kkt_matrix.Qxu_full().noalias() -= data.MJtJinv_dIDCdqv().transpose() * data.Qafu_full();
  kkt_residual.lx().noalias() -= data.MJtJinv_dIDCdqv().transpose() * data.laf();
  if (has_floating_base_) {
    kkt_matrix.Quu_full().noalias() = data.MJtJinv().topRows(dimv_) * data.Qafu_full();
    kkt_residual.lu_passive.noalias() += data.MJtJinv().template topRows<kDimFloatingBase>() * data.laf();
    kkt_residual.lu().noalias() += data.MJtJinv().block(kDimFloatingBase, 0, dimu_, dimv_+dimf_) * data.laf();
    kkt_residual.lu().noalias() -= kkt_matrix.Quu_drive_passive() * data.u_passive;
    kkt_residual.lx().noalias() -= kkt_matrix.Qxu_passive() * data.u_passive;
  }
  else {
    kkt_matrix.Quu().noalias() += data.MJtJinv().topRows(dimv_) * data.Qafu();
    kkt_residual.lu().noalias() += data.MJtJinv().topRows(dimv_) * data.laf();
  }
  kkt_matrix.Fvq() = - dtau * data.MJtJinv_dIDCdqv().topLeftCorner(dimv_, dimv_);
  kkt_matrix.Fvv().noalias() = Eigen::MatrixXd::Identity(dimv_, dimv_) 
                                - dtau * data.MJtJinv_dIDCdqv().topRightCorner(dimv_, dimv_);
  kkt_matrix.Fvu().noalias() += dtau * data.MJtJinv().block(0, dim_passive_, dimv_, dimu_);
  kkt_residual.Fv().noalias() -= dtau * data.MJtJinv_IDC().head(dimv_);
  if (has_floating_base_) {
    kkt_residual.Fv().noalias() -= dtau * data.MJtJinv().topLeftCorner(dimv_, dim_passive_) * data.u_passive;
  }
}


template <typename VectorType>
inline void ContactDynamics::computeCondensedDirection(
    const double dtau, const KKTMatrix& kkt_matrix, 
    const KKTResidual& kkt_residual, const Eigen::MatrixBase<VectorType>& dgmm, 
    SplitDirection& d) {
  assert(dtau > 0);
  assert(dgmm.size() == dimv_);
  if (has_floating_base_) {
    d.du_passive = - data_.u_passive;
    d.dnu_passive = kkt_residual.lu_passive;
    d.dnu_passive.noalias() += kkt_matrix.Quu_passive() * d.du_passive;
    d.dnu_passive.noalias() += kkt_matrix.Quu_passive_drive() * d.du();
    d.dnu_passive.noalias() += kkt_matrix.Qxu_passive().transpose() * d.dx();
    d.dnu_passive.noalias() 
        += dtau * data_.MJtJinv().topLeftCorner(dim_passive_, dimv_) * dgmm;
    d.dnu_passive.array() /= - dtau;
    d.daf().noalias() = - data_.MJtJinv_dIDCdqv() * d.dx();
    d.daf().noalias() += data_.MJtJinv().template leftCols<kDimFloatingBase>() * d.du_passive;
    d.daf().noalias() += data_.MJtJinv().block(0, kDimFloatingBase, dimv_+dimf_, dimu_) * d.du();
    d.daf().noalias() -= data_.MJtJinv_IDC();
  }
  else {
    d.daf().noalias() = - data_.MJtJinv_dIDCdqv() * d.dx();
    d.daf().noalias() += data_.MJtJinv().leftCols(dimv_) * d.du();
    d.daf().noalias() -= data_.MJtJinv_IDC();
  }
  d.df().array() *= -1;
  data_.laf().noalias() += data_.Qafqv() * d.dx();
  if (has_floating_base_) {
    data_.laf().noalias() += data_.Qafu().template leftCols<kDimFloatingBase>() * d.du_passive;
    data_.laf().noalias() += data_.Qafu().rightCols(dimv_-kDimFloatingBase) * d.du();
  }
  else {
    data_.laf().noalias() += data_.Qafu() * d.du();
  }
  data_.la().noalias() += dtau * dgmm;
  d.dbetamu().noalias() = - data_.MJtJinv() * data_.laf() / dtau;
}


inline void ContactDynamics::computeContactDynamicsResidual(
    Robot& robot, const ContactStatus& contact_status, const double dtau,
    const SplitSolution& s, KKTResidual& kkt_residual) {
  setContactStatus(contact_status);
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, ID_());
  if (has_floating_base_) {
    data_.ID().template head<kDimFloatingBase>().noalias() -= s.u_passive;
    data_.ID().tail(dimu_).noalias() -= s.u;
    u_passive_ = s.u_passive;
  }
  else {
    data_.ID().noalias() -= s.u;
  }
  robot.computeBaumgarteResidual(contact_status, dtau, data_.C());
}


inline void ContactDynamics::linearizeInverseDynamics(
    Robot& robot, const ContactStatus& contact_status, 
    const SplitSolution& s, KKTResidual& kkt_residual) {
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, data_.ID());
  if (has_floating_base_) {
    data_.ID().template head<kDimFloatingBase>().noalias() -= s.u_passive;
    data_.ID().tail(dimu_).noalias() -= s.u;
  }
  else {
    data_.ID().noalias() -= s.u;
  }
  robot.RNEADerivatives(s.q, s.v, s.a, data_.dIDdq(), data_.dIDdv(), data_.dIDda);
}


inline void ContactDynamics::linearizeContactConstraint(
    Robot& robot, const ContactStatus& contact_status, const double dtau,
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  assert(dtau > 0);
  robot.computeBaumgarteResidual(contact_status, dtau, data_.C());
  robot.computeBaumgarteDerivatives(contact_status, dtau, data_.dCdq(), 
                                    data_.dCdv(), data_.dCda());
}


inline double ContactDynamics::l1NormContactDynamicsResidual(
    const double dtau) const {
  assert(dtau > 0);
  return (dtau * (data_.IDC().lpNorm<1>() + data_.u_passive.lpNorm<1>()));
}


inline double ContactDynamics::squaredNormContactDynamicsResidual(
    const double dtau) const {
  assert(dtau > 0);
  return (dtau * dtau * (data_.IDC_().squaredNorm() + data_.u_passive.squaredNorm()));
}


inline void ContactDynamics::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
  has_active_contacts_ = contact_status.hasActiveContacts();
  data_.setContactStatus(contact_status);
}

} // namespace idocp 

#endif // IDOCP_CONTACT_DYNAMICS_HXX_ 