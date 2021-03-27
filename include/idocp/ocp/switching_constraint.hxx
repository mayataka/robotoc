#ifndef IDOCP_SWITCHING_CONSTRAINT_HXX_ 
#define IDOCP_SWITCHING_CONSTRAINT_HXX_

#include "idocp/ocp/switching_constraint.hpp"

#include <cassert>

namespace idocp {
namespace switchingconstraint {

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


inline double l1NormSwitchingConstraintResidual(
    const SplitKKTResidual& kkt_residual) {
  if (kkt_residual.P().size() > 0) { return kkt_residual.P().lpNorm<1>(); }
  else { return 0; }
}


inline double squaredNormSwitchingConstraintResidual(
    const SplitKKTResidual& kkt_residual) {
  if (kkt_residual.P().size() > 0) { return kkt_residual.P().squaredNorm(); }
  else { return 0; }
}

} // namespace switchingconstraint
} // namespace idocp

#endif // IDOCP_SWITCHING_CONSTRAINT_HXX_