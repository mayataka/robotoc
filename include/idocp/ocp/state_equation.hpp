#ifndef IDOCP_STATE_EQUATION_HPP_
#define IDOCP_STATE_EQUATION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {
namespace stateequation {

template <typename ConfigVectorType, typename SplitSolutionType>
void LinearizeForwardEuler(
    const Robot& robot, const double dtau, 
    const Eigen::MatrixBase<ConfigVectorType>& q_prev, const SplitSolution& s, 
    const SplitSolutionType& s_next, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual);

template <typename ConfigVectorType, typename TangentVectorType, 
          typename SplitSolutionType>
void LinearizeBackwardEuler(
    const Robot& robot, const double dtau, 
    const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const SplitSolution& s, const SplitSolutionType& s_next, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual);

template <typename ConfigVectorType, typename TangentVectorType>
void LinearizeBackwardEulerTerminal(
    const Robot& robot, const double dtau, 
    const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual);

template <typename ConfigVectorType, typename TangentVectorType>
void ComputeForwardEulerResidual(
    const Robot& robot, const double dtau, const SplitSolution& s, 
    const Eigen::MatrixBase<ConfigVectorType>& q_next, 
    const Eigen::MatrixBase<TangentVectorType>& v_next, 
    SplitKKTResidual& kkt_residual);

template <typename ConfigVectorType, typename TangentVectorType>
void ComputeBackwardEulerResidual(
    const Robot& robot, const double dtau, 
    const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual);

double L1NormStateEuqationResidual(const SplitKKTResidual& kkt_residual);

double SquaredNormStateEuqationResidual(const SplitKKTResidual& kkt_residual);

} // namespace stateequation

} // namespace idocp 

#include "idocp/ocp/state_equation.hxx"

#endif // IDOCP_STATE_EQUATION_HPP_