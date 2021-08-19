#ifndef IDOCP_SWITCHING_CONSTRAINT_HXX_ 
#define IDOCP_SWITCHING_CONSTRAINT_HXX_

#include "idocp/ocp/switching_constraint.hpp"

#include <cassert>

namespace idocp {
namespace switchingconstraint {

inline void linearizeSwitchingConstraint(
    Robot& robot, const ImpulseStatus& impulse_status, const double dt1, 
    const double dt2, const SplitSolution& s, SplitKKTResidual& kkt_residual, 
    SplitSwitchingConstraintJacobian& sc_jacobian,
    SplitSwitchingConstraintResidual& sc_residual) {
  assert(dt1 > 0);
  assert(dt2 > 0);
  sc_residual.setImpulseStatus(impulse_status);
  sc_jacobian.setImpulseStatus(impulse_status);
  sc_residual.setZero();
  sc_jacobian.setZero();
  computeSwitchingConstraintResidual(robot, impulse_status, dt1, dt2, s, 
                                     sc_residual);
  robot.computeContactPositionDerivative(impulse_status, sc_jacobian.Pq());
  if (robot.hasFloatingBase()) {
    robot.dIntegrateTransport_dq(s.q, sc_residual.dq, sc_jacobian.Pq(), 
                                 sc_jacobian.Phiq());
    robot.dIntegrateTransport_dv(s.q, sc_residual.dq, sc_jacobian.Pq(), 
                                 sc_jacobian.Phiv());
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
}


inline void computeSwitchingConstraintResidual(
    Robot& robot, const ImpulseStatus& impulse_status, const double dt1, 
    const double dt2, const SplitSolution& s, 
    SplitSwitchingConstraintResidual& sc_residual) {
  assert(dt1 > 0);
  assert(dt2 > 0);
  sc_residual.setImpulseStatus(impulse_status);
  sc_residual.dq = (dt1+dt2) * s.v + (dt1*dt2) * s.a;
  robot.integrateConfiguration(s.q, sc_residual.dq, 1.0, sc_residual.q);
  robot.updateKinematics(sc_residual.q);
  robot.computeContactPositionResidual(impulse_status, 
                                       impulse_status.contactPoints(), 
                                       sc_residual.P());
}

} // namespace switchingconstraint 
} // namespace idocp

#endif // IDOCP_SWITCHING_CONSTRAINT_HXX_