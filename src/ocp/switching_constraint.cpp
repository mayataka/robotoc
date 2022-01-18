#include "robotoc/ocp/switching_constraint.hpp"

#include <cassert>

namespace robotoc {

SwitchingConstraint::SwitchingConstraint(const Robot& robot)
  : q_(Eigen::VectorXd::Zero(robot.dimq())),
    dq_(Eigen::VectorXd::Zero(robot.dimv())),
    PqT_xi_(Eigen::VectorXd::Zero(robot.dimv())),
    has_floating_base_(robot.hasFloatingBase()) {
}


SwitchingConstraint::SwitchingConstraint()
  : q_(),
    dq_(),
    PqT_xi_(),
    has_floating_base_(false) {
}


SwitchingConstraint::~SwitchingConstraint() {
}


void SwitchingConstraint::evalSwitchingConstraint(
    Robot& robot, const ImpulseStatus& impulse_status, const double dt1, 
    const double dt2, const SplitSolution& s, 
    SwitchingConstraintResidual& sc_residual) {
  assert(dt1 > 0);
  assert(dt2 > 0);
  sc_residual.setImpulseStatus(impulse_status);
  dq_ = (dt1+dt2) * s.v + (dt1*dt2) * s.a;
  robot.integrateConfiguration(s.q, dq_, 1.0, q_);
  robot.updateKinematics(q_);
  robot.computeContactPositionResidual(impulse_status, sc_residual.P());
}


void SwitchingConstraint::linearizeSwitchingConstraint(
    Robot& robot, const ImpulseStatus& impulse_status, const double dt1, 
    const double dt2, const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual, SwitchingConstraintJacobian& sc_jacobian,
    SwitchingConstraintResidual& sc_residual) {
  assert(dt1 > 0);
  assert(dt2 > 0);
  sc_residual.setImpulseStatus(impulse_status);
  sc_jacobian.setImpulseStatus(impulse_status);
  sc_residual.setZero();
  sc_jacobian.setZero();
  evalSwitchingConstraint(robot, impulse_status, dt1, dt2, s, sc_residual);
  robot.computeContactPositionDerivative(impulse_status, sc_jacobian.Pq());
  if (has_floating_base_) {
    robot.dIntegrateTransport_dq(s.q, dq_, sc_jacobian.Pq(), sc_jacobian.Phiq());
    robot.dIntegrateTransport_dv(s.q, dq_, sc_jacobian.Pq(), sc_jacobian.Phiv());
    sc_jacobian.Phia() = (dt1*dt2) * sc_jacobian.Phiv();
    sc_jacobian.Phiv().array() *= (dt1+dt2);
  }
  else {
    sc_jacobian.Phiq() = sc_jacobian.Pq();
    sc_jacobian.Phiv() = (dt1+dt2) * sc_jacobian.Pq();
    sc_jacobian.Phia() = (dt1*dt2) * sc_jacobian.Pq();
  }
  kkt_residual.lx.noalias() += sc_jacobian.Phix().transpose() * s.xi_stack();
  kkt_residual.la.noalias() += sc_jacobian.Phia().transpose() * s.xi_stack();
  // STO sensitivities
  // Note that in computing the STO sensitivities, we always assume that dt1 = dt2.
  dq_ = 2.0 * (s.v + dt1 * s.a);
  sc_jacobian.Phit().noalias() = sc_jacobian.Pq() * dq_;
  kkt_residual.h += s.xi_stack().dot(sc_jacobian.Phit());
  PqT_xi_.noalias() = sc_jacobian.Pq().transpose() * s.xi_stack();
  kkt_matrix.Qtt += 2.0 * PqT_xi_.dot(s.a);
  kkt_matrix.hv().noalias() += 2.0 * PqT_xi_;
  kkt_matrix.ha.noalias() += (2.0*dt1) * PqT_xi_;
}

} // namespace robotoc