#ifndef IDOCP_IMPULSE_STATE_EQUATION_HPP_
#define IDOCP_IMPULSE_STATE_EQUATION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/state_equation.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"


namespace idocp {
namespace stateequation {

///
/// @brief Linearizes the impulse state equation of forward Euler. 
/// @param[in] robot Robot model. 
/// @param[in] q_prev Configuration at the previous time stage. 
/// @param[in] s Solution at the current impulse stage. 
/// @param[in] s_next Solution at the next time stage. 
/// @param[in, out] kkt_matrix Impulse split KKT matrix at the current impulse 
/// stage. 
/// @param[in, out] kkt_residual Impulse split KKT residual at the current 
/// impulse stage. 
///
template <typename ConfigVectorType>
void linearizeImpulseForwardEuler(
    const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const ImpulseSplitSolution& s, const SplitSolution& s_next, 
    ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual);

///
/// @brief Condenses terms related to the derivatives of the Lie group from the 
/// linearized impulse state equation of forward Euler. 
/// @param[in] robot Robot model. 
/// @param[in] s Solution at the current impulse stage. 
/// @param[in] q_next Configuration at the next time stage. 
/// @param[in, out] kkt_matrix Impulse split KKT matrix at the current impulse 
/// stage. 
/// @param[in, out] kkt_residual Impulse split KKT residual at the current 
/// impulse stage. 
///
template <typename ConfigVectorType>
void condenseImpulseForwardEuler(
    Robot& robot, const ImpulseSplitSolution& s, 
    const Eigen::MatrixBase<ConfigVectorType>& q_next, 
    ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual);

///
/// @brief Computes the residual in the impulse state equation of forward Euler. 
/// @param[in] robot Robot model. 
/// @param[in] s Solution at the current impulse stage. 
/// @param[in] q_next Configuration at the next time stage. 
/// @param[in] v_next Generalized velocity at the next time stage. 
/// @param[in, out] kkt_residual Impulse split KKT residual at the current 
/// impulse stage. 
///
template <typename ConfigVectorType, typename TangentVectorType>
void computeImpulseForwardEulerResidual(
    const Robot& robot, const ImpulseSplitSolution& s, 
    const Eigen::MatrixBase<ConfigVectorType>& q_next, 
    const Eigen::MatrixBase<TangentVectorType>& v_next, 
    ImpulseSplitKKTResidual& kkt_residual);

///
/// @brief Returns the l1-norm of the residual in the impulse state equation.
/// @param[in] kkt_residual Impulse split KKT residual at the current impulse stage. 
///
double l1NormStateEuqationResidual(
    const ImpulseSplitKKTResidual& kkt_residual);

///
/// @brief Returns the squared norm of the residual in the impulse state equation.
/// @param[in] kkt_residual Impulse split KKT residual at the current impulse stage. 
///
double squaredNormStateEuqationResidual(
    const ImpulseSplitKKTResidual& kkt_residual);

} // namespace stateequation
} // namespace idocp 

#include "idocp/impulse/impulse_state_equation.hxx"

#endif // IDOCP_IMPULSE_STATE_EQUATION_HPP_ 