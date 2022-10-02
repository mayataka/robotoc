#include "robotoc/dynamics/switching_constraint.hpp"

#include <cassert>

namespace robotoc {

void evalSwitchingConstraint(Robot& robot, const ImpulseStatus& impulse_status, 
                             SwitchingConstraintData& data, 
                             const double dt1, const double dt2, 
                             const SplitSolution& s, 
                             SplitKKTResidual& kkt_residual) {
  assert(dt1 > 0);
  assert(dt2 > 0);
  if (impulse_status.dimf() == 0) return;

  kkt_residual.setSwitchingConstraintDimension(impulse_status.dimf());
  kkt_residual.P().setZero();
  data.dq = (dt1+dt2) * s.v + (dt1*dt2) * s.a;
  robot.integrateConfiguration(s.q, data.dq, 1.0, data.q);
  robot.updateKinematics(data.q);
  robot.computeContactPositionResidual(impulse_status, kkt_residual.P());
}


void linearizeSwitchingConstraint(Robot& robot, 
                                  const ImpulseStatus& impulse_status, 
                                  SwitchingConstraintData& data, 
                                  const double dt1, const double dt2, 
                                  const SplitSolution& s, 
                                  SplitKKTMatrix& kkt_matrix, 
                                  SplitKKTResidual& kkt_residual) {
  assert(dt1 > 0);
  assert(dt2 > 0);
  if (impulse_status.dimf() == 0) return;

  kkt_matrix.setSwitchingConstraintDimension(impulse_status.dimf());
  kkt_residual.setSwitchingConstraintDimension(impulse_status.dimf());
  evalSwitchingConstraint(robot, impulse_status, data, dt1, dt2, s, kkt_residual);
  data.setDimension(impulse_status.dimf());
  data.Pq().setZero();
  robot.computeContactPositionDerivative(impulse_status, data.Pq());
  if (robot.hasFloatingBase()) {
    robot.dIntegrateTransport_dq(s.q, data.dq, data.Pq(), kkt_matrix.Phiq());
    robot.dIntegrateTransport_dv(s.q, data.dq, data.Pq(), kkt_matrix.Phiv());
    kkt_matrix.Phia() = (dt1*dt2) * kkt_matrix.Phiv();
    kkt_matrix.Phiv().array() *= (dt1+dt2);
  }
  else {
    kkt_matrix.Phiq() = data.Pq();
    kkt_matrix.Phiv() = (dt1+dt2) * data.Pq();
    kkt_matrix.Phia() = (dt1*dt2) * data.Pq();
  }
  kkt_residual.lx.noalias() += kkt_matrix.Phix().transpose() * s.xi_stack();
  kkt_residual.la.noalias() += kkt_matrix.Phia().transpose() * s.xi_stack();
  // STO sensitivities
  // Note that in computing the STO sensitivities, we always assume that dt1 = dt2.
  data.dq = 2.0 * (s.v + dt1 * s.a);
  kkt_matrix.Phit().noalias() = data.Pq() * data.dq;
  kkt_residual.h += s.xi_stack().dot(kkt_matrix.Phit());
  data.PqT_xi.noalias() = data.Pq().transpose() * s.xi_stack();
  kkt_matrix.Qtt += 2.0 * data.PqT_xi.dot(s.a);
  kkt_matrix.hv().noalias() += 2.0 * data.PqT_xi;
  kkt_matrix.ha.noalias() += (2.0*dt1) * data.PqT_xi;
}

} // namespace robotoc