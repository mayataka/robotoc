#ifndef IDOCP_CONTACT_DYNAMICS_HXX_
#define IDOCP_CONTACT_DYNAMICS_HXX_

#include "idocp/ocp/contact_dynamics.hpp"

#include <cassert>

namespace idocp {

inline ContactDynamics::ContactDynamics(const Robot& robot) 
  : data_(robot),
    has_floating_base_(robot.hasFloatingBase()),
    has_active_contacts_(false),
    dimv_(robot.dimv()),
    dimu_(robot.dimu()),
    dim_passive_(robot.dim_passive()) {
}


inline ContactDynamics::ContactDynamics() 
  : data_(),
    has_floating_base_(false),
    has_active_contacts_(false),
    dimv_(0),
    dimu_(0),
    dim_passive_(0) {
}


inline ContactDynamics::~ContactDynamics() {
}


inline void ContactDynamics::computeContactDynamicsResidual(
    Robot& robot, const ContactStatus& contact_status, const SplitSolution& s) {
  data_.setContactStatus(contact_status);
  has_active_contacts_ = contact_status.hasActiveContacts();
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, data_.ID_full());
  data_.ID().noalias() -= s.u;
  robot.computeBaumgarteResidual(contact_status, contact_status.contactPoints(), 
                                 data_.C());
}


inline void ContactDynamics::linearizeContactDynamics(
    Robot& robot, const ContactStatus& contact_status, const double dt, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) { 
  assert(dt > 0);
  computeContactDynamicsResidual(robot, contact_status, s);
  robot.RNEADerivatives(s.q, s.v, s.a, data_.dIDdq(), data_.dIDdv(), data_.dIDda);
  robot.computeBaumgarteDerivatives(contact_status, data_.dCdq(), data_.dCdv(), 
                                    data_.dCda());
  // augment inverse dynamics constraint
  kkt_residual.lq().noalias() += dt * data_.dIDdq().transpose() * s.beta;
  kkt_residual.lv().noalias() += dt * data_.dIDdv().transpose() * s.beta;
  kkt_residual.la.noalias()   += dt * data_.dIDda.transpose() * s.beta;
  if (has_active_contacts_) {
    kkt_residual.lf().noalias() -= dt * data_.dCda() * s.beta;
  }
  if (has_floating_base_) {
    // augment floating base constraint
    data_.lu_passive            = dt * s.nu_passive;
    data_.lu_passive.noalias() -= dt * s.beta.template head<kDimFloatingBase>(); 
    kkt_residual.lu.noalias()  -= dt * s.beta.tail(robot.dimu()); 
  }
  else {
    kkt_residual.lu.noalias() -= dt * s.beta; 
  }
  // augment acceleration-level contact constraint
  if (has_active_contacts_) {
    kkt_residual.lq().noalias() 
        += dt * data_.dCdq().transpose() * s.mu_stack();
    kkt_residual.lv().noalias() 
        += dt * data_.dCdv().transpose() * s.mu_stack();
    kkt_residual.la.noalias() += dt * data_.dCda().transpose() * s.mu_stack();
  }
}


inline void ContactDynamics::condenseContactDynamics(
    Robot& robot, const ContactStatus& contact_status, const double dt,
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dim_passive = robot.dim_passive();
  const int dimf = contact_status.dimf();
  robot.computeMJtJinv(data_.dIDda, data_.dCda(), data_.MJtJinv());
  data_.MJtJinv_dIDCdqv().noalias() = data_.MJtJinv() * data_.dIDCdqv();
  data_.MJtJinv_IDC().noalias()     = data_.MJtJinv() * data_.IDC();

  data_.Qafqv().topRows(dimv).noalias() 
      = (- kkt_matrix.Qaa.diagonal()).asDiagonal() 
          * data_.MJtJinv_dIDCdqv().topRows(dimv);
  data_.Qafqv().bottomRows(dimf).noalias() 
      = - kkt_matrix.Qff() * data_.MJtJinv_dIDCdqv().bottomRows(dimf);
  data_.Qafqv().bottomLeftCorner(dimf, dimv).noalias()
      -= kkt_matrix.Qqf().transpose();
  data_.Qafu_full().topRows(dimv).noalias() 
      = kkt_matrix.Qaa.diagonal().asDiagonal() 
          * data_.MJtJinv().topLeftCorner(dimv, dimv);
  data_.Qafu_full().bottomRows(dimf).noalias() 
      = kkt_matrix.Qff() * data_.MJtJinv().bottomLeftCorner(dimf, dimv);
  data_.la() = kkt_residual.la;
  data_.lf() = - kkt_residual.lf();
  data_.la().noalias() 
      -= kkt_matrix.Qaa.diagonal().asDiagonal() 
          * data_.MJtJinv_IDC().head(dimv);
  data_.lf().noalias() 
      -= kkt_matrix.Qff() * data_.MJtJinv_IDC().tail(dimf);

  kkt_matrix.Qxx.noalias() 
      -= data_.MJtJinv_dIDCdqv().transpose() * data_.Qafqv();
  kkt_matrix.Qxx.topRows(dimv).noalias() 
      += kkt_matrix.Qqf() * data_.MJtJinv_dIDCdqv().bottomRows(dimf);
  if (has_floating_base_) {
    data_.Qxu_passive.noalias() 
        = - data_.MJtJinv_dIDCdqv().transpose() * data_.Qafu_full().leftCols(dim_passive);
    data_.Qxu_passive.topRows(dimv).noalias()
        -= kkt_matrix.Qqf() * data_.MJtJinv().bottomLeftCorner(dimf, dimv).leftCols(dim_passive);
    kkt_matrix.Qxu.noalias() 
        -= data_.MJtJinv_dIDCdqv().transpose() * data_.Qafu_full().rightCols(dimu);
    kkt_matrix.Qxu.topRows(dimv).noalias()
        -= kkt_matrix.Qqf() * data_.MJtJinv().bottomLeftCorner(dimf, dimv).rightCols(dimu);
  }
  else {
    kkt_matrix.Qxu.noalias() 
        -= data_.MJtJinv_dIDCdqv().transpose() * data_.Qafu_full();
    kkt_matrix.Qxu.topRows(dimv).noalias()
        -= kkt_matrix.Qqf() * data_.MJtJinv().bottomLeftCorner(dimf, dimv);
  }
  kkt_residual.lx.noalias() 
      -= data_.MJtJinv_dIDCdqv().transpose() * data_.laf();
  kkt_residual.lq().noalias()
      += kkt_matrix.Qqf() * data_.MJtJinv_IDC().tail(dimf);

  if (has_floating_base_) {
    data_.Quu_passive_topRight.noalias() 
        = data_.MJtJinv().topRows(dim_passive) * data_.Qafu_full().rightCols(dimu);
    kkt_matrix.Quu.noalias() 
        += data_.MJtJinv().middleRows(dim_passive, dimu) * data_.Qafu_full().rightCols(dimu);
  }
  else {
    kkt_matrix.Quu.noalias() 
        += data_.MJtJinv().topRows(dimv) * data_.Qafu_full();
  }
  if (has_floating_base_) {
    data_.lu_passive.noalias() 
        += data_.MJtJinv().template topRows<kDimFloatingBase>() * data_.laf();
  }
  kkt_residual.lu.noalias() 
      += data_.MJtJinv().middleRows(dim_passive, dimu) * data_.laf();

  kkt_matrix.Fvq() = - dt * data_.MJtJinv_dIDCdqv().topLeftCorner(dimv, dimv);
  kkt_matrix.Fvv().noalias() 
        = - dt * data_.MJtJinv_dIDCdqv().topRightCorner(dimv, dimv) 
          + Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix.Fvu = dt * data_.MJtJinv().block(0, dim_passive, dimv, dimu);
  kkt_residual.Fv().noalias() -= dt * data_.MJtJinv_IDC().head(dimv);
}


inline void ContactDynamics::expandPrimal(SplitDirection& d) const {
  d.daf().noalias() = - data_.MJtJinv_dIDCdqv() * d.dx;
  d.daf().noalias() 
      += data_.MJtJinv().middleCols(dim_passive_, dimu_) * d.du;
  d.daf().noalias() -= data_.MJtJinv_IDC();
  d.df().array()    *= -1;
}


template <typename SplitDirectionType>
inline void ContactDynamics::expandDual(const double dt, 
                                        const SplitDirectionType& d_next, 
                                        SplitDirection& d) {
  assert(dt > 0);
  if (has_floating_base_) {
    d.dnu_passive            = data_.lu_passive;
    d.dnu_passive.noalias() += data_.Quu_passive_topRight * d.du;
    d.dnu_passive.noalias() += data_.Qxu_passive.transpose() * d.dx;
    d.dnu_passive.noalias() 
        += dt * data_.MJtJinv().leftCols(dimv_).template topRows<kDimFloatingBase>() 
              * d_next.dgmm();
    d.dnu_passive.array() *= - (1/dt);
  }
  data_.laf().noalias() += data_.Qafqv() * d.dx;
  data_.laf().noalias() += data_.Qafu() * d.du;
  data_.la().noalias()  += dt * d_next.dgmm();
  d.dbetamu().noalias()  = - data_.MJtJinv() * data_.laf() * (1/dt);
}


inline void ContactDynamics::condenseSwitchingConstraint(
    SplitSwitchingConstraintJacobian& sc_jacobian,
    SplitSwitchingConstraintResidual& sc_residual) const {
  sc_jacobian.Phix().noalias() 
      -= sc_jacobian.Phia() * data_.MJtJinv_dIDCdqv().topRows(dimv_);
  sc_jacobian.Phiu().noalias()  
      = sc_jacobian.Phia() * data_.MJtJinv().block(0, dim_passive_, dimv_, dimu_);
  sc_residual.P().noalias() 
      -= sc_jacobian.Phia() * data_.MJtJinv_IDC().topRows(dimv_);
}


inline double ContactDynamics::KKTError() const {
  return (data_.IDC().squaredNorm() + data_.lu_passive.squaredNorm());
}


inline double ContactDynamics::constraintViolation() const {
  return data_.IDC().lpNorm<1>();
}

} // namespace idocp 

#endif // IDOCP_CONTACT_DYNAMICS_HXX_ 