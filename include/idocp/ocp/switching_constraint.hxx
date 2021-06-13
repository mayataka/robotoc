#ifndef IDOCP_SWITCHING_CONSTRAINT_HXX_ 
#define IDOCP_SWITCHING_CONSTRAINT_HXX_

#include "idocp/ocp/switching_constraint.hpp"

#include <cassert>

namespace idocp {
namespace switchingconstraint {

inline void linearizeSwitchingConstraint(
    Robot& robot, const ImpulseStatus& impulse_status, const double dt1, 
    const double dt2, const SplitSolution& s, SplitKKTResidual& kkt_residual, 
    SplitSwitchingConstraintJacobian& switching_jacobian,
    SplitSwitchingConstraintResidual& switching_residual) {
  assert(dt1 > 0);
  assert(dt2 > 0);
  switching_residual.setImpulseStatus(impulse_status);
  switching_jacobian.setImpulseStatus(impulse_status);
  switching_residual.setZero();
  switching_jacobian.setZero();
  computeSwitchingConstraintResidual(robot, impulse_status, dt1, dt2, s, 
                                     switching_residual);
  robot.computeContactPositionDerivative(impulse_status, switching_jacobian.Pq());
  if (robot.hasFloatingBase()) {
    robot.dIntegrateTransport_dq(s.q, switching_residual.dq, 
                                 switching_jacobian.Pq(), 
                                 switching_jacobian.Phiq());
    robot.dIntegrateTransport_dv(s.q, switching_residual.dq, 
                                 switching_jacobian.Pq(), 
                                 switching_jacobian.Phiv());
    switching_jacobian.Phia() = (dt1*dt2) * switching_jacobian.Phiv();
    switching_jacobian.Phiv().array() *= (dt1+dt2);
  }
  else {
    switching_jacobian.Phiq() = switching_jacobian.Pq();
    switching_jacobian.Phiv() = (dt1+dt2) * switching_jacobian.Pq();
    switching_jacobian.Phia() = (dt1*dt2) * switching_jacobian.Pq();
  }
  kkt_residual.lx.noalias() 
      += switching_jacobian.Phix().transpose() * s.xi_stack();
  kkt_residual.la.noalias() 
      += switching_jacobian.Phia().transpose() * s.xi_stack();
}


inline void computeSwitchingConstraintResidual(
    Robot& robot, const ImpulseStatus& impulse_status, const double dt1, 
    const double dt2, const SplitSolution& s, 
    SplitSwitchingConstraintResidual& switching_residual) {
  assert(dt1 > 0);
  assert(dt2 > 0);
  switching_residual.setImpulseStatus(impulse_status);
  switching_residual.dq = (dt1+dt2) * s.v + (dt1*dt2) * s.a;
  robot.integrateConfiguration(s.q, switching_residual.dq, 1.0, 
                               switching_residual.q);
  robot.updateKinematics(switching_residual.q);
  robot.computeContactPositionResidual(impulse_status, 
                                       impulse_status.contactPoints(), 
                                       switching_residual.P());
}

} // namespace switchingconstraint 
} // namespace idocp

#endif // IDOCP_SWITCHING_CONSTRAINT_HXX_