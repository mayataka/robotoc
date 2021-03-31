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

///
/// @brief Linearizes the switching constraint, i.e., the contact position
/// constraint according to the input impulse status, for the backward Euler.
/// @param[in] robot Robot model. Kinematics must be updated.
/// @param[in] impulse_status Impulse status. 
/// @param[in] s Split solution of the time stage just before the impulse.
/// @param[in, out] kkt_matrix Split KKT matrix of the time stage just before 
/// the impulse.
/// @param[in, out] kkt_residual Split KKT residual of the time stage just 
/// before the impulse.
///
void linearizeSwitchingConstraint(Robot& robot, 
                                  const ImpulseStatus& impulse_status, 
                                  const SplitSolution& s, 
                                  SplitKKTMatrix& kkt_matrix, 
                                  SplitKKTResidual& kkt_residual);

///
/// @brief Computes the residual in the switching constraint, i.e., the contact 
/// position constraint according to the input impulse status, for the backward 
/// Euler.
/// @param[in] robot Robot model. Kinematics must be updated.
/// @param[in] impulse_status Impulse status. 
/// @param[in, out] kkt_residual Split KKT residual of the time stage just 
/// before the impulse.
///
void computeSwitchingConstraintResidual(Robot& robot, 
                                        const ImpulseStatus& impulse_status,
                                        SplitKKTResidual& kkt_residual);

///
/// @brief Returns l1-norm of the residual in the switching constraint. 
/// @param[in] kkt_residual Split KKT residual of the time stage just before
/// the impulse.
/// @return l1-norm of the residual in the switching constraint.
///
double l1NormSwitchingConstraintResidual(const SplitKKTResidual& kkt_residual);

///
/// @brief Returns squared norm of the residual in the switching constraint. 
/// @param[in] kkt_residual Split KKT residual of the time stage just before
/// the impulse.
/// @return Squared norm of the residual in the switching constraint.
///
double squaredNormSwitchingConstraintResidual(
    const SplitKKTResidual& kkt_residual);

} // namespace switchingconstraint
} // namespace idocp

#include "idocp/ocp/switching_constraint.hxx"

#endif // IDOCP_SWITCHING_CONSTRAINT_HPP_