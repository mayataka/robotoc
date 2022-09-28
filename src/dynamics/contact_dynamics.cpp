#include "robotoc/dynamics/contact_dynamics.hpp"

#include <cassert>


namespace robotoc {

namespace {
  constexpr int dim_floating_base = 6;
} 

void evalContactDynamics(Robot& robot, const ContactStatus& contact_status, 
                         const SplitSolution& s, ContactDynamicsData& data) {
  data.setContactDimension(contact_status.dimf());
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, data.ID_full());
  data.ID().noalias() -= s.u;
  robot.computeBaumgarteResidual(contact_status, data.C());
}


void linearizeContactDynamics(Robot& robot, const ContactStatus& contact_status, 
                              const SplitSolution& s, 
                              ContactDynamicsData& data, 
                              SplitKKTResidual& kkt_residual) { 
  evalContactDynamics(robot, contact_status, s, data);
  robot.RNEADerivatives(s.q, s.v, s.a, data.dIDdq(), data.dIDdv(), data.dIDda);
  robot.computeBaumgarteDerivatives(contact_status, data.dCdq(), data.dCdv(), 
                                    data.dCda());
  // augment inverse dynamics constraint
  kkt_residual.lq().noalias() += data.dIDdq().transpose() * s.beta;
  kkt_residual.lv().noalias() += data.dIDdv().transpose() * s.beta;
  kkt_residual.la.noalias()   += data.dIDda.transpose() * s.beta;
  if (contact_status.hasActiveContacts()) {
    kkt_residual.lf().noalias() -= data.dCda() * s.beta;
  }
  if (data.hasFloatingBase()) {
    // augment floating base constraint
    data.lu_passive            = s.nu_passive;
    data.lu_passive.noalias() -= s.beta.template head<dim_floating_base>(); 
    kkt_residual.lu.noalias()  -= s.beta.tail(robot.dimu()); 
  }
  else {
    kkt_residual.lu.noalias() -= s.beta; 
  }
  // augment acceleration-level contact constraint
  if (contact_status.hasActiveContacts()) {
    kkt_residual.lq().noalias() += data.dCdq().transpose() * s.mu_stack();
    kkt_residual.lv().noalias() += data.dCdv().transpose() * s.mu_stack();
    kkt_residual.la.noalias()   += data.dCda().transpose() * s.mu_stack();
  }
}


void condenseContactDynamics(Robot& robot, const ContactStatus& contact_status, 
                             const double dt, ContactDynamicsData& data, 
                             SplitKKTMatrix& kkt_matrix, 
                             SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dim_passive = robot.dim_passive();
  const int dimf = contact_status.dimf();
  robot.computeMJtJinv(data.dIDda, data.dCda(), data.MJtJinv());
  data.MJtJinv_dIDCdqv().noalias() = data.MJtJinv() * data.dIDCdqv();
  data.MJtJinv_IDC().noalias()     = data.MJtJinv() * data.IDC();

  data.Qafqv().topRows(dimv).noalias() 
      = (- kkt_matrix.Qaa.diagonal()).asDiagonal() 
          * data.MJtJinv_dIDCdqv().topRows(dimv);
  data.Qafqv().bottomRows(dimf).noalias() 
      = - kkt_matrix.Qff() * data.MJtJinv_dIDCdqv().bottomRows(dimf);
  data.Qafqv().bottomLeftCorner(dimf, dimv).noalias()
      -= kkt_matrix.Qqf().transpose();
  data.Qafu_full().topRows(dimv).noalias() 
      = kkt_matrix.Qaa.diagonal().asDiagonal() 
          * data.MJtJinv().topLeftCorner(dimv, dimv);
  data.Qafu_full().bottomRows(dimf).noalias() 
      = kkt_matrix.Qff() * data.MJtJinv().bottomLeftCorner(dimf, dimv);
  data.la() = kkt_residual.la;
  data.lf() = - kkt_residual.lf();
  data.la().noalias() 
      -= kkt_matrix.Qaa.diagonal().asDiagonal() 
          * data.MJtJinv_IDC().head(dimv);
  data.lf().noalias() 
      -= kkt_matrix.Qff() * data.MJtJinv_IDC().tail(dimf);

  kkt_matrix.Qxx.noalias() 
      -= data.MJtJinv_dIDCdqv().transpose() * data.Qafqv();
  kkt_matrix.Qxx.topRows(dimv).noalias() 
      += kkt_matrix.Qqf() * data.MJtJinv_dIDCdqv().bottomRows(dimf);
  if (data.hasFloatingBase()) {
    data.Qxu_passive.noalias() 
        = - data.MJtJinv_dIDCdqv().transpose() * data.Qafu_full().leftCols(dim_passive);
    data.Qxu_passive.topRows(dimv).noalias()
        -= kkt_matrix.Qqf() * data.MJtJinv().bottomLeftCorner(dimf, dimv).leftCols(dim_passive);
    kkt_matrix.Qxu.noalias() 
        -= data.MJtJinv_dIDCdqv().transpose() * data.Qafu_full().rightCols(dimu);
    kkt_matrix.Qxu.topRows(dimv).noalias()
        -= kkt_matrix.Qqf() * data.MJtJinv().bottomLeftCorner(dimf, dimv).rightCols(dimu);
  }
  else {
    kkt_matrix.Qxu.noalias() 
        -= data.MJtJinv_dIDCdqv().transpose() * data.Qafu_full();
    kkt_matrix.Qxu.topRows(dimv).noalias()
        -= kkt_matrix.Qqf() * data.MJtJinv().bottomLeftCorner(dimf, dimv);
  }
  kkt_residual.lx.noalias() 
      -= data.MJtJinv_dIDCdqv().transpose() * data.laf();
  kkt_residual.lq().noalias()
      += kkt_matrix.Qqf() * data.MJtJinv_IDC().tail(dimf);

  if (data.hasFloatingBase()) {
    data.Quu_passive_topRight.noalias() 
        = data.MJtJinv().topRows(dim_passive) * data.Qafu_full().rightCols(dimu);
    kkt_matrix.Quu.noalias() 
        += data.MJtJinv().middleRows(dim_passive, dimu) * data.Qafu_full().rightCols(dimu);
  }
  else {
    kkt_matrix.Quu.noalias() 
        += data.MJtJinv().topRows(dimv) * data.Qafu_full();
  }
  if (data.hasFloatingBase()) {
    data.lu_passive.noalias() 
        += data.MJtJinv().template topRows<dim_floating_base>() * data.laf();
  }
  kkt_residual.lu.noalias() 
      += data.MJtJinv().middleRows(dim_passive, dimu) * data.laf();

  kkt_matrix.Fvq() = - dt * data.MJtJinv_dIDCdqv().topLeftCorner(dimv, dimv);
  kkt_matrix.Fvv().noalias() 
        = - dt * data.MJtJinv_dIDCdqv().topRightCorner(dimv, dimv) 
          + Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix.Fvu = dt * data.MJtJinv().block(0, dim_passive, dimv, dimu);
  kkt_residual.Fv().noalias() -= dt * data.MJtJinv_IDC().head(dimv);

  // Switching constraint
  if (kkt_matrix.dims() > 0) {
    assert(kkt_matrix.dims() == kkt_residual.dims());
    data.setSwitchingConstraintDimension(kkt_matrix.dims());
    data.Phia() = kkt_matrix.Phia();
    kkt_matrix.Phix().noalias() 
        -= data.Phia() * data.MJtJinv_dIDCdqv().topRows(data.dimv());
    kkt_matrix.Phiu().noalias()  
        = data.Phia() * data.MJtJinv().block(0, data.dim_passive(), data.dimv(), data.dimu());
    kkt_matrix.Phit().noalias() 
        -= data.Phia() * data.MJtJinv_IDC().head(data.dimv());
    kkt_residual.P().noalias() 
        -= data.Phia() * data.MJtJinv_IDC().head(data.dimv());
  }
  else {
    data.setSwitchingConstraintDimension(0);
  }

  // STO sensitivities
  data.ha() = kkt_matrix.ha;
  data.hf() = - kkt_matrix.hf();
  kkt_residual.h -= data.MJtJinv_IDC().dot(data.haf()); 
  kkt_matrix.hx.noalias() -= data.MJtJinv_dIDCdqv().transpose() * data.haf();
  kkt_matrix.hq().noalias()
      += (1.0/dt) * kkt_matrix.Qqf() * data.MJtJinv_IDC().tail(dimf);
  kkt_matrix.hu.noalias() 
      += data.MJtJinv().middleRows(dim_passive, dimu) * data.haf();
}


void condenseContactDynamics(const ContactDynamicsData& data,
                             SwitchingConstraintJacobian& sc_jacobian,
                             SwitchingConstraintResidual& sc_residual) {
  sc_jacobian.Phix().noalias() 
      -= sc_jacobian.Phia() * data.MJtJinv_dIDCdqv().topRows(data.dimv());
  sc_jacobian.Phiu().noalias()  
      = sc_jacobian.Phia() * data.MJtJinv().block(0, data.dim_passive(), data.dimv(), data.dimu());
  sc_jacobian.Phit().noalias() 
      -= sc_jacobian.Phia() * data.MJtJinv_IDC().head(data.dimv());
  sc_residual.P().noalias() 
      -= sc_jacobian.Phia() * data.MJtJinv_IDC().head(data.dimv());
}


void expandContactDynamicsPrimal(const ContactDynamicsData& data, 
                                 SplitDirection& d) {
  d.daf().noalias() = - data.MJtJinv_dIDCdqv() * d.dx;
  d.daf().noalias() 
      += data.MJtJinv().middleCols(data.dim_passive(), data.dimu()) * d.du;
  d.daf().noalias() -= data.MJtJinv_IDC();
  d.df().array()    *= -1;
}


void expandContactDynamicsDual(const double dt, const double dts, 
                               ContactDynamicsData& data, 
                               const SplitDirection& d_next, 
                               SplitDirection& d) {
  assert(dt > 0);
  if (data.hasFloatingBase()) {
    d.dnu_passive            = - data.lu_passive;
    d.dnu_passive.noalias() -= data.Quu_passive_topRight * d.du;
    d.dnu_passive.noalias() -= data.Qxu_passive.transpose() * d.dx;
    d.dnu_passive.noalias() 
        -= dt * data.MJtJinv().leftCols(data.dimv()).template topRows<dim_floating_base>() 
              * d_next.dgmm();
  }
  data.laf().noalias() += data.Qafqv() * d.dx;
  data.laf().noalias() += data.Qafu() * d.du;
  data.la().noalias()  += dt * d_next.dgmm();
  if (data.dims() > 0) {
    assert(data.dims() == d.dims());
    data.la().noalias()  += data.Phia().transpose() * d.dxi();
  }
  constexpr double eps = std::numeric_limits<double>::epsilon();
  if (dts < - eps || dts > eps) {
    data.laf().noalias() += dts * data.haf();
  }
  d.dbetamu().noalias()  = - data.MJtJinv() * data.laf();
}


void expandContactDynamicsDual(const double dt, const double dts, 
                               ContactDynamicsData& data, 
                               const SwitchingConstraintJacobian& sc_jacobian,
                               const SplitDirection& d_next, SplitDirection& d) {
  assert(dt > 0);
  if (data.hasFloatingBase()) {
    d.dnu_passive            = - data.lu_passive;
    d.dnu_passive.noalias() -= data.Quu_passive_topRight * d.du;
    d.dnu_passive.noalias() -= data.Qxu_passive.transpose() * d.dx;
    d.dnu_passive.noalias() 
        -= dt * data.MJtJinv().leftCols(data.dimv()).template topRows<dim_floating_base>() 
              * d_next.dgmm();
  }
  data.laf().noalias() += data.Qafqv() * d.dx;
  data.laf().noalias() += data.Qafu() * d.du;
  data.la().noalias()  += dt * d_next.dgmm();
  data.la().noalias()  += sc_jacobian.Phia().transpose() * d.dxi();
  constexpr double eps = std::numeric_limits<double>::epsilon();
  if (dts < - eps || dts > eps) {
    data.laf().noalias() += dts * data.haf();
  }
  d.dbetamu().noalias()  = - data.MJtJinv() * data.laf();
}

} // namespace robotoc 