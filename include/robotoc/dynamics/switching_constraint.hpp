#ifndef ROBOTOC_SWITCHING_CONSTRAINT_HPP_ 
#define ROBOTOC_SWITCHING_CONSTRAINT_HPP_

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impact_status.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/dynamics/switching_constraint_data.hpp"


namespace robotoc {

///
/// @brief Computes the residual in the switching constraint, i.e., the 
/// contact position constraint. 
/// @note The internal kinematics data of robot is updated.  
/// @param[in] robot Robot model. Kinematics must be updated.
/// @param[in] impact_status Impact status. 
/// @param[in, out] data Data structure for the switching constraint. 
/// @param[in] dt1 Time step of the time stage 2 stage before the impact.
/// @param[in] dt2 Time step of the time stage just before the impact.
/// @param[in] s Split solution of the time stage 2 stage before the impact.
/// @param[in, out] kkt_residual Split KKT residual of the time stage 2 stage
/// before the impact.
///
void evalSwitchingConstraint(Robot& robot, const ImpactStatus& impact_status, 
                             SwitchingConstraintData& data, 
                             const double dt1, const double dt2, 
                             const SplitSolution& s, 
                             SplitKKTResidual& kkt_residual);

///
/// @brief Linearizes the switching constraint, i.e., the contact position 
/// constraint. 
/// @note The internal kinematics data of robot is updated.  
/// @param[in] robot Robot model. Kinematics must be updated.
/// @param[in] impact_status Impact status. 
/// @param[in, out] data Data structure for the switching constraint. 
/// @param[in] dt1 Time step of the time stage 2 stage before the impact.
/// @param[in] dt2 Time step of the time stage just before the impact.
/// @param[in] s Split solution of the time stage 2 stage before the impact.
/// @param[in, out] kkt_matrix Split KKT matrix of the time stage 2 stage
/// before the impact.
/// @param[in, out] kkt_residual Split KKT residual of the time stage 2 stage
/// before the impact.
///
void linearizeSwitchingConstraint(
    Robot& robot, const ImpactStatus& impact_status, 
    SwitchingConstraintData& data, const double dt1, const double dt2, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual);

} // namespace robotoc

#endif // ROBOTOC_SWITCHING_CONSTRAINT_HPP_ 