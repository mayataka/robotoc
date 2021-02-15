#ifndef IDOCP_CONTACT_DYNAMICS_HXX_
#define IDOCP_CONTACT_DYNAMICS_HXX_

#include "idocp/ocp/contact_dynamics.hpp"

#include <stdexcept>
#include <cassert>

namespace idocp {

inline ContactDynamics::ContactDynamics(const Robot& robot, 
                                        const double baumgarte_time_step) 
  : data_(robot),
    has_floating_base_(robot.hasFloatingBase()),
    has_active_contacts_(false),
    baumgarte_time_step_(baumgarte_time_step) {
  try {
    if (baumgarte_time_step <= 0) {
      throw std::out_of_range(  
          "invalid value: baumgarte_time_step must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline ContactDynamics::ContactDynamics() 
  : data_(),
    has_floating_base_(false),
    has_active_contacts_(false),
    baumgarte_time_step_(0) {
}


inline ContactDynamics::~ContactDynamics() {
}


inline void ContactDynamics::linearizeContactDynamics(
    Robot& robot, const ContactStatus& contact_status, const double dtau, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) { 
  assert(dtau >= 0);
  setContactStatus(contact_status);
  linearizeInverseDynamics(robot, contact_status, s, data_);
  linearizeContactConstraint(robot, contact_status, baumgarte_time_step_, data_);
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
    kkt_residual.lu_passive = dtau * s.nu_passive;
    kkt_residual.lu_passive.noalias() -= dtau * s.beta.template head<kDimFloatingBase>(); 
    kkt_residual.lu().noalias() -= dtau * s.beta.tail(robot.dimu()); 
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
  data.ID().noalias() -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, data.dIDdq(), data.dIDdv(), data.dIDda);
}


inline void ContactDynamics::linearizeContactConstraint(
    Robot& robot, const ContactStatus& contact_status, 
    const double baumgarte_time_step, ContactDynamicsData& data) {
  assert(baumgarte_time_step > 0);
  robot.computeBaumgarteResidual(contact_status, baumgarte_time_step, 
                                 contact_status.contactPoints(), data.C());
  robot.computeBaumgarteDerivatives(contact_status, baumgarte_time_step, 
                                    data.dCdq(), data.dCdv(), data.dCda());
}


inline void ContactDynamics::condenseContactDynamics(
    Robot& robot, const ContactStatus& contact_status, const double dtau,
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual,
    const bool is_forward_euler) {
  assert(dtau >= 0);
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dim_passive = robot.dim_passive();
  const int dimf = kkt_matrix.dimf();
  robot.computeMJtJinv(data_.dIDda, data_.dCda(), data_.MJtJinv());
  data_.MJtJinv_dIDCdqv().noalias() = data_.MJtJinv() * data_.dIDCdqv();
  data_.MJtJinv_IDC().noalias() = data_.MJtJinv() * data_.IDC();
  data_.Qafqv().topRows(dimv).noalias() 
      = (- kkt_matrix.Qaa().diagonal()).asDiagonal() 
          * data_.MJtJinv_dIDCdqv().topRows(dimv);
  data_.Qafqv().bottomRows(dimf).noalias() 
      = - kkt_matrix.Qff() * data_.MJtJinv_dIDCdqv().bottomRows(dimf);
  data_.Qafu_full().topRows(dimv).noalias() 
      = kkt_matrix.Qaa().diagonal().asDiagonal() 
          * data_.MJtJinv().leftCols(dimv).topRows(dimv);
  data_.Qafu_full().bottomRows(dimf).noalias() 
      = kkt_matrix.Qff() * data_.MJtJinv().leftCols(dimv).bottomRows(dimf);
  data_.la() = kkt_residual.la;
  data_.lf() = - kkt_residual.lf();
  data_.la().noalias() 
      -= kkt_matrix.Qaa().diagonal().asDiagonal() 
          * data_.MJtJinv_IDC().head(dimv);
  data_.lf().noalias() 
      -= kkt_matrix.Qff() * data_.MJtJinv_IDC().tail(dimf);
  kkt_matrix.Qxx().noalias() 
      -= data_.MJtJinv_dIDCdqv().transpose() * data_.Qafqv();
  kkt_matrix.Qxu_full().noalias() 
      -= data_.MJtJinv_dIDCdqv().transpose() * data_.Qafu_full();
  kkt_residual.lx().noalias() 
      -= data_.MJtJinv_dIDCdqv().transpose() * data_.laf();
  kkt_matrix.Quu_full().noalias() 
      += data_.MJtJinv().topRows(dimv) * data_.Qafu_full();
  if (has_floating_base_) {
    kkt_residual.lu_passive.noalias() 
        += data_.MJtJinv().template topRows<kDimFloatingBase>() * data_.laf();
  }
  kkt_residual.lu().noalias() 
      += data_.MJtJinv().middleRows(dim_passive, dimu) * data_.laf();
  kkt_matrix.Fvq() = - dtau * data_.MJtJinv_dIDCdqv().topLeftCorner(dimv, dimv);
  if (is_forward_euler) {
    kkt_matrix.Fvv().noalias() 
        = - dtau * data_.MJtJinv_dIDCdqv().topRightCorner(dimv, dimv) 
          + Eigen::MatrixXd::Identity(dimv, dimv);
  }
  else {
    kkt_matrix.Fvv().noalias() 
        = - dtau * data_.MJtJinv_dIDCdqv().topRightCorner(dimv, dimv) 
          - Eigen::MatrixXd::Identity(dimv, dimv);
  }
  kkt_matrix.Fvu() = dtau * data_.MJtJinv().block(0, dim_passive, dimv, dimu);
  kkt_residual.Fv().noalias() -= dtau * data_.MJtJinv_IDC().head(dimv);
}


inline void ContactDynamics::computeCondensedPrimalDirection(
    const Robot& robot, SplitDirection& d) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  d.daf().noalias() = - data_.MJtJinv_dIDCdqv() * d.dx();
  d.daf().noalias() 
      += data_.MJtJinv().middleCols(robot.dim_passive(), dimu) * d.du();
  d.daf().noalias() -= data_.MJtJinv_IDC();
  d.df().array() *= -1;
}


template <typename VectorType>
inline void ContactDynamics::computeCondensedDualDirection(
    const Robot& robot, const double dtau, const SplitKKTMatrix& kkt_matrix, 
    const SplitKKTResidual& kkt_residual, 
    const Eigen::MatrixBase<VectorType>& dgmm, SplitDirection& d) {
  assert(dtau >= 0);
  assert(dgmm.size() == robot.dimv());
  if (dtau >= kMindtau) {
    const int dimv = robot.dimv();
    if (has_floating_base_) {
      d.dnu_passive = kkt_residual.lu_passive;
      d.dnu_passive.noalias() += kkt_matrix.Quu_passive_topRight() * d.du();
      d.dnu_passive.noalias() += kkt_matrix.Qxu_passive().transpose() * d.dx();
      d.dnu_passive.noalias() 
          += dtau * data_.MJtJinv().leftCols(dimv).template topRows<kDimFloatingBase>() * dgmm;
      d.dnu_passive.array() *= - (1/dtau);
    }
    data_.laf().noalias() += data_.Qafqv() * d.dx();
    data_.laf().noalias() += data_.Qafu() * d.du();
    data_.la().noalias() += dtau * dgmm;
    d.dbetamu().noalias() = - data_.MJtJinv() * data_.laf() * (1/dtau);
  }
  else {
    // In this case (dtau < kMindtau) we regard dtau as zero. Then we set the 
    // directions of the dual variables zero because they are undefined.
    if (has_floating_base_) {
      d.dnu_passive.setZero();
    }
    d.dbetamu().setZero();
  }
}


inline void ContactDynamics::computeContactDynamicsResidual(
    Robot& robot, const ContactStatus& contact_status, const SplitSolution& s) {
  setContactStatus(contact_status);
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, data_.ID_full());
  data_.ID().noalias() -= s.u;
  robot.computeBaumgarteResidual(contact_status, baumgarte_time_step_, 
                                 contact_status.contactPoints(), data_.C());
}


inline double ContactDynamics::l1NormContactDynamicsResidual(
    const double dtau) const {
  assert(dtau >= 0);
  return (dtau * data_.IDC().lpNorm<1>());
}


inline double ContactDynamics::squaredNormContactDynamicsResidual(
    const double dtau) const {
  assert(dtau >= 0);
  return (dtau * dtau * data_.IDC().squaredNorm());
}


inline void ContactDynamics::setContactStatus(
    const ContactStatus& contact_status) {
  data_.setContactStatus(contact_status);
  has_active_contacts_ = contact_status.hasActiveContacts();
}

} // namespace idocp 

#endif // IDOCP_CONTACT_DYNAMICS_HXX_ 