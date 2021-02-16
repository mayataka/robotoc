#ifndef IDOCP_SWITCHING_CONSTRAINT_HPP_ 
#define IDOCP_SWITCHING_CONSTRAINT_HPP_

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_state_constraint_jacobian.hpp"


namespace idocp {
namespace switchingconstraint {

void linearizeSwitchingConstraint(Robot& robot, 
                                  const ImpulseStatus& impulse_status, 
                                  const double dtau1, const double dtau2, 
                                  const SplitSolution& s, 
                                  SplitKKTMatrix& kkt_matrix, 
                                  SplitKKTResidual& kkt_residual, 
                                  SplitStateConstraintJacobian& jac);

void computeSwitchingConstraintResidual(Robot& robot, 
                                        const ImpulseStatus& impulse_status, 
                                        const double dtau1, const double dtau2, 
                                        const SplitSolution& s, 
                                        SplitKKTResidual& kkt_residual, 
                                        SplitStateConstraintJacobian& jac);

void linearizeSwitchingConstraint(Robot& robot, 
                                  const ImpulseStatus& impulse_status, 
                                  const SplitSolution& s, 
                                  SplitKKTMatrix& kkt_matrix, 
                                  SplitKKTResidual& kkt_residual);

void computeSwitchingConstraintResidual(Robot& robot, 
                                        const ImpulseStatus& impulse_status,
                                        SplitKKTResidual& kkt_residual);

double squaredNormSwitchingConstraintResidual(
    const SplitKKTResidual& kkt_residual);

double l1NormSwitchingConstraintResidual(const SplitKKTResidual& kkt_residual);

} // namespace switchingconstraint
} // namespace idocp

#include "idocp/ocp/switching_constraint.hxx"

#endif // IDOCP_SWITCHING_CONSTRAINT_HPP_