#include "robotoc/dynamics/switching_constraint.hpp"

#include <cassert>

namespace robotoc {

SwitchingConstraint::SwitchingConstraint(const Robot& robot)
  : data_(robot),
    has_floating_base_(robot.hasFloatingBase()) {
}


SwitchingConstraint::SwitchingConstraint()
  : data_(),
    has_floating_base_(false) {
}


void SwitchingConstraint::evalSwitchingConstraint(
    Robot& robot, const ImpulseStatus& impulse_status, const double dt1, 
    const double dt2, const SplitSolution& s, 
    SwitchingConstraintResidual& sc_residual) {
  assert(dt1 > 0);
  assert(dt2 > 0);
  sc_residual.setDimension(impulse_status.dimf());
  data_.dq = (dt1+dt2) * s.v + (dt1*dt2) * s.a;
  robot.integrateConfiguration(s.q, data_.dq, 1.0, data_.q);
  robot.updateKinematics(data_.q);
  robot.computeContactPositionResidual(impulse_status, sc_residual.P());
}


void SwitchingConstraint::linearizeSwitchingConstraint(
    Robot& robot, const ImpulseStatus& impulse_status, const double dt1, 
    const double dt2, const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual, SwitchingConstraintJacobian& sc_jacobian,
    SwitchingConstraintResidual& sc_residual) {
  assert(dt1 > 0);
  assert(dt2 > 0);
  sc_residual.setDimension(impulse_status.dimf());
  sc_jacobian.setDimension(impulse_status.dimf());
  sc_residual.setZero();
  sc_jacobian.setZero();
  evalSwitchingConstraint(robot, impulse_status, dt1, dt2, s, sc_residual);
  robot.computeContactPositionDerivative(impulse_status, sc_jacobian.Pq());
  if (has_floating_base_) {
    robot.dIntegrateTransport_dq(s.q, data_.dq, sc_jacobian.Pq(), sc_jacobian.Phiq());
    robot.dIntegrateTransport_dv(s.q, data_.dq, sc_jacobian.Pq(), sc_jacobian.Phiv());
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
  data_.dq = 2.0 * (s.v + dt1 * s.a);
  sc_jacobian.Phit().noalias() = sc_jacobian.Pq() * data_.dq;
  kkt_residual.h += s.xi_stack().dot(sc_jacobian.Phit());
  data_.PqT_xi.noalias() = sc_jacobian.Pq().transpose() * s.xi_stack();
  kkt_matrix.Qtt += 2.0 * data_.PqT_xi.dot(s.a);
  kkt_matrix.hv().noalias() += 2.0 * data_.PqT_xi;
  kkt_matrix.ha.noalias() += (2.0*dt1) * data_.PqT_xi;
}

} // namespace robotoc