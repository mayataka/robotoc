#ifndef IDOCP_SWITCHING_CONSTRAINT_HXX_ 
#define IDOCP_SWITCHING_CONSTRAINT_HXX_

#include "idocp/ocp/switching_constraint.hpp"

#include <cassert>

namespace idocp {
namespace switchingconstraint {

inline void linearizeSwitchingConstraint(
    Robot& robot, const ImpulseStatus& impulse_status, const double dt1, 
    const double dt2, const SplitSolution& s, SplitKKTResidual& kkt_residual, 
    SplitSwitchingConstraintJacobian& switch_jacobian,
    SplitSwitchingConstraintResidual& switch_residual) {
  assert(dt1 > 0);
  assert(dt2 > 0);
  switch_residual.setImpulseStatus(impulse_status);
  switch_jacobian.setImpulseStatus(impulse_status);
  switch_residual.setZero();
  switch_jacobian.setZero();
  computeSwitchingConstraintResidual(robot, impulse_status, dt1, dt2, s, 
                                     switch_residual);
  robot.computeContactPositionDerivative(impulse_status, switch_jacobian.Pq());
  if (robot.hasFloatingBase()) {
    robot.dIntegrateTransport_dq(s.q, switch_residual.dq_pred, 
                                 switch_jacobian.Pq(), switch_jacobian.Phiq());
    robot.dIntegrateTransport_dv(s.q, switch_residual.dq_pred, 
                                 switch_jacobian.Pq(), switch_jacobian.Phiv());
    switch_jacobian.Phia() = (dt1*dt2) * switch_jacobian.Phiv();
    switch_jacobian.Phiv().array() *= (dt1+dt2);
  }
  else {
    switch_jacobian.Phiq() = switch_jacobian.Pq();
    switch_jacobian.Phiv() = (dt1+dt2) * switch_jacobian.Pq();
    switch_jacobian.Phia() = (dt1*dt2) * switch_jacobian.Pq();
  }
  kkt_residual.lx.noalias() += switch_jacobian.Phix().transpose() * s.xi_stack();
  kkt_residual.la.noalias() += switch_jacobian.Phia().transpose() * s.xi_stack();
}


inline void computeSwitchingConstraintResidual(
    Robot& robot, const ImpulseStatus& impulse_status, const double dt1, 
    const double dt2, const SplitSolution& s, 
    SplitSwitchingConstraintResidual& switch_residual) {
  assert(dt1 > 0);
  assert(dt2 > 0);
  switch_residual.setImpulseStatus(impulse_status);
  switch_residual.dq_pred = (dt1+dt2) * s.v + (dt1*dt2) * s.a;
  robot.integrateConfiguration(s.q, switch_residual.dq_pred, 1.0, 
                               switch_residual.q_pred);
  robot.updateKinematics(switch_residual.q_pred);
  robot.computeContactPositionResidual(impulse_status, 
                                       impulse_status.contactPoints(), 
                                       switch_residual.P());
}


inline double l1NormSwitchingConstraintResidual(
    const SplitSwitchingConstraintResidual& switch_residual) {
  return switch_residual.P().template lpNorm<1>();
}


inline double squaredNormSwitchingConstraintResidual(
    const SplitSwitchingConstraintResidual& switch_residual) {
  return switch_residual.P().squaredNorm();
}

} // namespace switchingconstraint 
} // namespace idocp

#endif // IDOCP_SWITCHING_CONSTRAINT_HXX_