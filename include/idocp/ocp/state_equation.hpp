#ifndef IDOCP_STATE_EQUATION_HPP_
#define IDOCP_STATE_EQUATION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {
namespace stateequation {

///
/// @brief Linearizes the state equation of forward Euler. 
/// @param[in] robot Robot model. 
/// @param[in] dt Time step. 
/// @param[in] q_prev Configuration at the previous time stage. 
/// @param[in] s Solution at the current stage. 
/// @param[in] s_next Solution at the next time stage. 
/// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
/// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
///
template <typename ConfigVectorType, typename SplitSolutionType>
void linearizeForwardEuler(
    const Robot& robot, const double dt, 
    const Eigen::MatrixBase<ConfigVectorType>& q_prev, const SplitSolution& s, 
    const SplitSolutionType& s_next, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual);


///
/// @brief Condenses terms related to the derivatives of the Lie group from the 
/// linearized state equation of forward Euler. 
/// @param[in] robot Robot model. 
/// @param[in] dt Time step. 
/// @param[in] s Solution at the current time stage. 
/// @param[in] q_next Configuration at the next time stage. 
/// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
/// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
///
template <typename ConfigVectorType>
void condenseForwardEuler(
    Robot& robot, const double dt, const SplitSolution& s, 
    const Eigen::MatrixBase<ConfigVectorType>& q_next, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual);

///
/// @brief Linearizes the state equation of forward Euler at the teminal stage. 
/// @param[in] robot Robot model. 
/// @param[in] q_prev Configuration at the previous time stage. 
/// @param[in] s Solution at the current stage. 
/// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
/// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
///
template <typename ConfigVectorType>
void linearizeForwardEulerTerminal(
    const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual);

///
/// @brief Condenses terms related to the derivatives of the Lie group from the 
/// linearized state equation of forward Euler at the terminal stage. 
/// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
///
void condenseForwardEulerTerminal(Robot& robot, SplitKKTMatrix& kkt_matrix);

///
/// @brief Corrects the costate direction using the derivatives of the Lie group. 
/// @param[in] robot Robot model. 
/// @param[in] kkt_matrix Split KKT matrix at the current time stage. 
/// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
/// @param[in, out] dlmd Costate direction. 
///
template <typename SplitKKTMatrixType, typename SplitKKTResidualType, 
          typename VectorType>
void correctCostateDirectionForwardEuler(
    const Robot& robot, const SplitKKTMatrixType& kkt_matrix, 
    SplitKKTResidualType& kkt_residual,
    const Eigen::MatrixBase<VectorType>& dlmd);

///
/// @brief Linearizes the state equation of backward Euler. 
/// @param[in] robot Robot model. 
/// @param[in] dt Time step. 
/// @param[in] q_prev Configuration at the previous time stage. 
/// @param[in] v_prev Generalized velocity at the previous time stage. 
/// @param[in] s Solution at the current tiem stage. 
/// @param[in] s_next Solution at the next time stage. 
/// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
/// @param[in, out] kkt_residual Split KKT reisdual at the current time stage. 
///
template <typename ConfigVectorType, typename TangentVectorType, 
        typename SplitSolutionType>
void linearizeBackwardEuler(
    const Robot& robot, const double dt, 
    const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const SplitSolution& s, const SplitSolutionType& s_next, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual);

///
/// @brief Linearizes the state equation of backward Euler at the terminal stage. 
/// @param[in] robot Robot model. 
/// @param[in] dt Time step. 
/// @param[in] q_prev Configuration at the previous time stage. 
/// @param[in] v_prev Generalized velocity at the previous time stage. 
/// @param[in] s Solution at the current tiem stage. 
/// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
/// @param[in, out] kkt_residual Split KKT reisdual at the current time stage. 
///
template <typename ConfigVectorType, typename TangentVectorType>
void linearizeBackwardEulerTerminal(
    const Robot& robot, const double dt, 
    const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual);

///
/// @brief Condenses terms related to the derivatives of the Lie group from the 
/// linearized state equation of backward Euler. 
/// @param[in] robot Robot model. 
/// @param[in] dt Time step. 
/// @param[in] q_prev Configuration at the previous time stage. 
/// @param[in] s Solution at the current tiem stage. 
/// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
/// @param[in, out] kkt_residual Split KKT reisdual at the current time stage. 
///
template <typename ConfigVectorType>
void condenseBackwardEuler(Robot& robot, const double dt, 
                           const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
                           const SplitSolution& s, 
                           SplitKKTMatrix& kkt_matrix, 
                           SplitKKTResidual& kkt_residual);

///
/// @brief Corrects the costate direction using the derivatives of the Lie group. 
/// @param[in] robot Robot model. 
/// @param[in] kkt_matrix Split KKT matrix at the current time stage. 
/// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
/// @param[in, out] dlmd Costate direction. 
///
template <typename SplitKKTMatrixType, typename SplitKKTResidualType, 
          typename VectorType>
void correctCostateDirectionBackwardEuler(
    const Robot& robot, const SplitKKTMatrixType& kkt_matrix, 
    SplitKKTResidualType& kkt_residual,
    const Eigen::MatrixBase<VectorType>& dlmd);

///
/// @brief Computes the residual in the state equation of forward Euler. 
/// @param[in] robot Robot model. 
/// @param[in] dt Time step. 
/// @param[in] s Solution at the current time stage. 
/// @param[in] q_next Configuration at the next time stage. 
/// @param[in] v_next Generalized velocity at the next time stage. 
/// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
///
template <typename ConfigVectorType, typename TangentVectorType>
void computeForwardEulerResidual(
    const Robot& robot, const double dt, const SplitSolution& s, 
    const Eigen::MatrixBase<ConfigVectorType>& q_next, 
    const Eigen::MatrixBase<TangentVectorType>& v_next, 
    SplitKKTResidual& kkt_residual);

///
/// @brief Computes the residual in the state equation of backward Euler. 
/// @param[in] robot Robot model. 
/// @param[in] dt Time step. 
/// @param[in] q_prev Configuration at the previous time stage. 
/// @param[in] v_prev Generalized velocity at the previous time stage. 
/// @param[in] s Solution at the current time stage. 
/// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
///
template <typename ConfigVectorType, typename TangentVectorType>
void computeBackwardEulerResidual(
    const Robot& robot, const double dt, 
    const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual);

///
/// @brief Returns the l1-norm of the residual in the state equation.
/// @param[in] kkt_residual Split KKT residual at the current time stage. 
///
double l1NormStateEuqationResidual(
    const SplitKKTResidual& kkt_residual);

///
/// @brief Returns the squared norm of the residual in the state equation.
/// @param[in] kkt_residual Split KKT residual at the current time stage. 
///
double squaredNormStateEuqationResidual(
    const SplitKKTResidual& kkt_residual);

} // namespace stateequation 
} // namespace idocp 

#include "idocp/ocp/state_equation.hxx"

#endif // IDOCP_STATE_EQUATION_HPP_