#ifndef ROBOTOC_CONSTRAINTS_IMPL_HPP_
#define ROBOTOC_CONSTRAINTS_IMPL_HPP_

#include <vector>
#include <memory>

#include "robotoc/robot/robot.hpp"
#include "robotoc/constraints/constraint_component_base.hpp"
#include "robotoc/constraints/impulse_constraint_component_base.hpp"
#include "robotoc/constraints/constraint_component_data.hpp"
#include "robotoc/constraints/constraints_data.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/impulse/impulse_split_solution.hpp"
#include "robotoc/impulse/impulse_split_direction.hpp"
#include "robotoc/impulse/impulse_split_kkt_residual.hpp"
#include "robotoc/impulse/impulse_split_kkt_matrix.hpp"


namespace robotoc {
namespace constraintsimpl {

using ConstraintComponentBasePtr = std::shared_ptr<ConstraintComponentBase>;
using ImpulseConstraintComponentBasePtr 
  = std::shared_ptr<ImpulseConstraintComponentBase>;

///
/// @brief Clears the vector of the constraints. 
/// @param[in, out] constraints Vector of the constraints. 
///
template <typename ConstraintComponentBaseTypePtr>
void clear(std::vector<ConstraintComponentBaseTypePtr>& constraints);

///
/// @brief Checks whether the constrainst use kinematics or not. 
/// @param[in] constraints Vector of the constraints. 
/// @return true if the constraints requre kinematics of Robot model. false if not.
///
template <typename ConstraintComponentBaseTypePtr>
bool useKinematics(
   const std::vector<ConstraintComponentBaseTypePtr>& constraints);

///
/// @brief Creates constraints data.
/// @param[in] constraints Vector of the constraints. 
/// @param[out] data Vector of the constraints data. 
///
template <typename ConstraintComponentBaseTypePtr>
void createConstraintsData(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints, 
    std::vector<ConstraintComponentData>& data);

///
/// @brief Checks whether the current solution s is feasible or not. 
/// @param[in] constraints Vector of the constraints. 
/// @param[in] robot Robot model.
/// @param[in, out] data Vector of the constraints data. 
/// @param[in] s Split solution.
/// @return true if s is feasible. false if not.
///
template <typename ConstraintComponentBaseTypePtr, typename SplitSolutionType>
bool isFeasible(const std::vector<ConstraintComponentBaseTypePtr>& constraints,
                Robot& robot, std::vector<ConstraintComponentData>& data, 
                const SplitSolutionType& s);

///
/// @brief Sets the slack and dual variables of each constraint components. 
/// @param[in] constraints Vector of the constraints. 
/// @param[in] robot Robot model.
/// @param[in, out] data Vector of the constraints data. 
/// @param[in] s Split solution.
///
template <typename ConstraintComponentBaseTypePtr, typename SplitSolutionType>
void setSlackAndDual(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const SplitSolutionType& s);

///
/// @brief Computes the primal residual, residual in the complementary 
/// slackness, and the log-barrier function of the slack varible.
/// @param[in] constraints Vector of the impulse constraints. 
/// @param[in] robot Robot model.
/// @param[in, out] data Vector of the constraints data.
/// @param[in] s Split solution.
///
template <typename ConstraintComponentBaseTypePtr, typename SplitSolutionType>
void evalConstraint(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const SplitSolutionType& s);

///
/// @brief Evaluates the constraints (i.e., calls evalConstraint()) and adds 
/// the products of the Jacobian of the constraints and Lagrange multipliers.
/// @param[in] constraints Vector of the impulse constraints. 
/// @param[in] robot Robot model.
/// @param[in, out] data Vector of the constraints data.
/// @param[in] dt Time step.
/// @param[in] s Impulse split solution.
/// @param[in, out] kkt_residual Impulse split KKT residual.
///
void evalDerivatives(
    const std::vector<ConstraintComponentBasePtr>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, const double dt,
    const SplitSolution& s, SplitKKTResidual& kkt_residual);

///
/// @brief Evaluates the constraints (i.e., calls evalConstraint()) and adds 
/// the products of the Jacobian of the constraints and Lagrange multipliers.
/// @param[in] constraints Vector of the impulse constraints. 
/// @param[in] robot Robot model.
/// @param[in, out] data Vector of the constraints data.
/// @param[in] s Impulse split solution.
/// @param[in, out] kkt_residual Impulse split KKT residual.
///
void evalDerivatives(
    const std::vector<ImpulseConstraintComponentBasePtr>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTResidual& kkt_residual);

///
/// @brief Condenses the slack and dual variables. linearizeConstraints() must 
/// be called before this function.
/// @param[in] constraints Vector of the constraints. 
/// @param[in] robot Robot model.
/// @param[in, out] data Vector of the constraints data.
/// @param[in] dt Time step.
/// @param[in] s Split solution.
/// @param[in, out] kkt_matrix Split KKT matrix. The condensed Hessians are added  
/// to this object.
/// @param[in, out] kkt_residual Split KKT residual. The condensed residuals are 
/// added to this object.
///
void condenseSlackAndDual(
    const std::vector<ConstraintComponentBasePtr>& constraints, Robot& robot, 
    std::vector<ConstraintComponentData>& data, const double dt, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual);

///
/// @brief Condenses the slack and dual variables. linearizeConstraints() must 
/// be called before this function.
/// @param[in] constraints Vector of the impulse constraints. 
/// @param[in] robot Robot model.
/// @param[in, out] data Vector of the constraints data.
/// @param[in] s Impulse split solution.
/// @param[in, out] kkt_matrix Impulse split KKT matrix. The condensed Hessians are 
/// added to this object.
/// @param[in, out] kkt_residual Impulse split KKT residual. The condensed residuals 
/// are added to this object.
///
void condenseSlackAndDual(
    const std::vector<ImpulseConstraintComponentBasePtr>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual);

///
/// @brief Expands the slack and dual, i.e., computes the directions of the 
/// slack and dual variables from the directions of the primal variables.
/// @param[in] constraints Vector of the impulse constraints. 
/// @param[in, out] data Vector of the constraints data.
/// @param[in] dt Time step of this time stage. 
/// @param[in] dts The direction of the switching time regarding of this time 
/// stage. 
/// @param[in] s Split solution.
/// @param[in] d Split direction.
///
template <typename ConstraintComponentBaseTypePtr, 
          typename SplitSolutionType, typename SplitDirectionType>
void expandSlackAndDual(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    std::vector<ConstraintComponentData>& data, const double dt, 
    const double dts, const SplitSolutionType& s, const SplitDirectionType& d);

///
/// @brief Expands the slack and dual, i.e., computes the directions of the 
/// slack and dual variables from the directions of the primal variables.
/// @param[in] constraints Vector of the impulse constraints. 
/// @param[in, out] data Vector of the constraints data.
/// @param[in] s Split solution.
/// @param[in] d Split direction.
///
template <typename ConstraintComponentBaseTypePtr, 
          typename SplitSolutionType, typename SplitDirectionType>
void expandSlackAndDual(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    std::vector<ConstraintComponentData>& data, 
    const SplitSolutionType& s, const SplitDirectionType& d);

///
/// @brief Computes and returns the maximum step size by applying 
/// fraction-to-boundary-rule to the direction of the slack variable.
/// @param[in] constraints Vector of the impulse constraints. 
/// @param[in] data Vector of the constraints data.
///
template <typename ConstraintComponentBaseTypePtr>
double maxSlackStepSize(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    const std::vector<ConstraintComponentData>& data);

///
/// @brief Computes and returns the maximum step size by applying 
/// fraction-to-boundary-rule to the direction of the dual variable.
/// @param[in] constraints Vector of the impulse constraints. 
/// @param[in] data Vector of the constraints data.
///
template <typename ConstraintComponentBaseTypePtr>
double maxDualStepSize(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    const std::vector<ConstraintComponentData>& data);

///
/// @brief Updates the slack variables according to the step size.
/// @param[in, out] data Vector of the constraints data.
/// @param[in] step_size Step size. 
///
void updateSlack(std::vector<ConstraintComponentData>& data, 
                 const double step_size);

///
/// @brief Updates the dual variables according to the step size.
/// @param[in, out] data Vector of the constraints data.
/// @param[in] step_size Step size. 
///
void updateDual(std::vector<ConstraintComponentData>& data, 
                const double step_size);

///
/// @brief Sets the barrier parameter.
/// @param[in, out] constraints Vector of the impulse constraints. 
/// @param[in] barrier Barrier parameter. Must be positive. Should be small.
///
template <typename ConstraintComponentBaseTypePtr>
void setBarrier(std::vector<ConstraintComponentBaseTypePtr>& constraints, 
                const double barrier);

///
/// @brief Sets the parameter of the fraction-to-boundary-rule.
/// @param[in, out] constraints Vector of the impulse constraints. 
/// @param[in] fraction_to_boundary_rule Must be larger than 0 and smaller 
/// than 1. Should be between 0.9 and 0.995.
///
template <typename ConstraintComponentBaseTypePtr>
void setFractionToBoundaryRule(
    std::vector<ConstraintComponentBaseTypePtr>& constraints,
    const double fraction_to_boundary_rule);

} // namespace constraintsimpl
} // namespace robotoc

#include "robotoc/constraints/constraints_impl.hxx"

#endif // ROBOTOC_CONSTRAINTS_IMPL_HPP_