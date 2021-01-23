#ifndef IDOCP_STATE_EQUATION_HPP_
#define IDOCP_STATE_EQUATION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {
namespace stateequation {
    
void condenseForwardEuler(Robot& robot, const double dtau, 
                          SplitKKTMatrix& kkt_matrix, 
                          SplitKKTResidual& kkt_residual);

void condenseForwardEulerTerminal(Robot& robot, SplitKKTMatrix& kkt_matrix);

template <typename TangentVectorType>
void correctCostateDirectionForwardEuler(
    const Robot& robot, const SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual,
    const Eigen::MatrixBase<TangentVectorType>& dlmd);

void condenseBackwardEuler(Robot& robot, const double dtau, 
                           SplitKKTMatrix& kkt_matrix, 
                           SplitKKTResidual& kkt_residual);

template <typename ConfigVectorType, typename SplitSolutionType>
void linearizeForwardEuler(
    const Robot& robot, const double dtau, 
    const Eigen::MatrixBase<ConfigVectorType>& q_prev, const SplitSolution& s, 
    const SplitSolutionType& s_next, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual);

template <typename ConfigVectorType>
void linearizeForwardEulerTerminal(
    const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual);

template <typename ConfigVectorType, typename TangentVectorType, 
        typename SplitSolutionType>
void linearizeBackwardEuler(
    const Robot& robot, const double dtau, 
    const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const SplitSolution& s, const SplitSolutionType& s_next, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual);

template <typename ConfigVectorType, typename TangentVectorType>
void linearizeBackwardEulerTerminal(
    const Robot& robot, const double dtau, 
    const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual);

template <typename ConfigVectorType, typename TangentVectorType>
void computeForwardEulerResidual(
    const Robot& robot, const double dtau, const SplitSolution& s, 
    const Eigen::MatrixBase<ConfigVectorType>& q_next, 
    const Eigen::MatrixBase<TangentVectorType>& v_next, 
    SplitKKTResidual& kkt_residual);

template <typename ConfigVectorType, typename TangentVectorType>
void computeBackwardEulerResidual(
    const Robot& robot, const double dtau, 
    const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual);

double l1NormStateEuqationResidual(
    const SplitKKTResidual& kkt_residual);

double squaredNormStateEuqationResidual(
    const SplitKKTResidual& kkt_residual);

} // namespace stateequation 
} // namespace idocp 

#include "idocp/ocp/state_equation.hxx"

#endif // IDOCP_STATE_EQUATION_HPP_