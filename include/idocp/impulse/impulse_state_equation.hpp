#ifndef IDOCP_IMPULSE_STATE_EQUATION_HPP_
#define IDOCP_IMPULSE_STATE_EQUATION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/state_equation.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"


namespace idocp {
namespace stateequation {

template <typename ConfigVectorType>
void LinearizeImpulseForwardEuler(
    const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const ImpulseSplitSolution& s, const SplitSolution& s_next, 
    ImpulseKKTMatrix& kkt_matrix, 
    ImpulseKKTResidual& kkt_residual);

template <typename ConfigVectorType, typename TangentVectorType>
void LinearizeImpulseBackwardEuler(
    const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const ImpulseSplitSolution& s, const SplitSolution& s_next, 
    ImpulseKKTMatrix& kkt_matrix, ImpulseKKTResidual& kkt_residual);

template <typename ConfigVectorType, typename TangentVectorType>
void LinearizeImpulseBackwardEulerTerminal(
    const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, KKTResidual& kkt_residual);

template <typename ConfigVectorType, typename TangentVectorType>
void ComputeImpulseForwardEulerResidual(
    const Robot& robot, const ImpulseSplitSolution& s, 
    const Eigen::MatrixBase<ConfigVectorType>& q_next, 
    const Eigen::MatrixBase<TangentVectorType>& v_next, 
    ImpulseKKTResidual& kkt_residual);

template <typename ConfigVectorType, typename TangentVectorType>
void ComputeImpulseBackwardEulerResidual(
    const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const ImpulseSplitSolution& s, ImpulseKKTResidual& kkt_residual);

double L1NormStateEuqationResidual(const ImpulseKKTResidual& kkt_residual);

double SquaredNormStateEuqationResidual(const ImpulseKKTResidual& kkt_residual);

} // namespace stateequation

} // namespace idocp 

#include "idocp/impulse/impulse_state_equation.hxx"

#endif // IDOCP_IMPULSE_STATE_EQUATION_HPP_ 