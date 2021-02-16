#ifndef IDOCP_SWITCHING_CONSTRAINT_HXX_ 
#define IDOCP_SWITCHING_CONSTRAINT_HXX_

#include "idocp/ocp/switching_constraint.hpp"

#include <cassert>

namespace idocp {
namespace switchingconstraint {

inline void linearizeSwitchingConstraint(
    Robot& robot, const ImpulseStatus& impulse_status, const double dtau1, 
    const double dtau2, const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual, SplitStateConstraintJacobian& jac) {
  assert(dtau1 > 0);
  assert(dtau2 > 0);
  jac.setImpulseStatus(impulse_status);
  computeSwitchingConstraintResidual(robot, impulse_status, dtau1, dtau2, s, 
                                     kkt_residual, jac);
  robot.computeContactDerivative(impulse_status, kkt_matrix.Pq());
  if (robot.hasFloatingBase()) {
    robot.dIntegratedConfiguration(s.q, jac.dq, jac.dintegrate_dq);
    robot.dIntegratedVelocity(s.q, jac.dq, jac.dintegrate_dv);
    jac.Phiq().noalias() = kkt_matrix.Pq() * jac.dintegrate_dq;
    jac.Phiv().noalias() = (dtau1+dtau2) * kkt_matrix.Pq() * jac.dintegrate_dv;
    jac.Phia().noalias() = (dtau1*dtau2) * kkt_matrix.Pq() * jac.dintegrate_dv;
  }
  else {
    jac.Phiq() = kkt_matrix.Pq();
    jac.Phiv() = (dtau1+dtau2) * kkt_matrix.Pq();
    jac.Phia() = (dtau1*dtau2) * kkt_matrix.Pq();
  }
  kkt_residual.lq().noalias() += jac.Phiq().transpose() * s.xi_stack();
  kkt_residual.lv().noalias() += jac.Phiv().transpose() * s.xi_stack();
  kkt_residual.la.noalias()   += jac.Phia().transpose() * s.xi_stack();
}


inline void computeSwitchingConstraintResidual(
    Robot& robot, const ImpulseStatus& impulse_status, const double dtau1, 
    const double dtau2, const SplitSolution& s, SplitKKTResidual& kkt_residual, 
    SplitStateConstraintJacobian& jac) {
  assert(dtau1 > 0);
  assert(dtau2 > 0);
  jac.dq = (dtau1+dtau2) * s.v + (dtau1*dtau2) * s.a;
  robot.integrateConfiguration(s.q, jac.dq, 1.0, jac.q);
  robot.updateKinematics(jac.q);
  robot.computeContactResidual(impulse_status, impulse_status.contactPoints(), 
                               kkt_residual.P());
}


inline void linearizeSwitchingConstraint(
    Robot& robot, const ImpulseStatus& impulse_status, const SplitSolution& s, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) {
  computeSwitchingConstraintResidual(robot, impulse_status, kkt_residual);
  robot.computeContactDerivative(impulse_status, kkt_matrix.Pq());
  kkt_residual.lq().noalias() += kkt_matrix.Pq().transpose() * s.xi_stack();
}


inline void computeSwitchingConstraintResidual(
    Robot& robot, const ImpulseStatus& impulse_status,
    SplitKKTResidual& kkt_residual) {
  robot.computeContactResidual(impulse_status, impulse_status.contactPoints(),
                               kkt_residual.P());
}


inline double squaredNormSwitchingConstraintResidual(
    const SplitKKTResidual& kkt_residual) {
  return kkt_residual.P().squaredNorm();
}


inline double l1NormSwitchingConstraintResidual(
    const SplitKKTResidual& kkt_residual) {
  return kkt_residual.P().lpNorm<1>();
}

} // namespace switchingconstraint
} // namespace idocp

#endif // IDOCP_SWITCHING_CONSTRAINT_HXX_