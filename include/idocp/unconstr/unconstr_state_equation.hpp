#ifndef IDOCP_UNCONSTR_STATE_EQUATION_HPP_
#define IDOCP_UNCONSTR_STATE_EQUATION_HPP_

#include "Eigen/Core"

#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {
namespace unconstr {
namespace stateequation {

///
/// @brief Linearizes the state equation of forward Euler. 
/// @param[in] dt Time step. 
/// @param[in] s Solution at the current stage. 
/// @param[in] s_next Solution at the next time stage. 
/// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
/// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
///
void linearizeForwardEuler(const double dt, const SplitSolution& s, 
                           const SplitSolution& s_next, 
                           SplitKKTMatrix& kkt_matrix, 
                           SplitKKTResidual& kkt_residual);

///
/// @brief Linearizes the state equation of backward Euler. 
/// @param[in] dt Time step. 
/// @param[in] q_prev Configuration at the previous time stage. 
/// @param[in] v_prev Generalized velocity at the previous time stage. 
/// @param[in] s Solution at the current tiem stage. 
/// @param[in] s_next Solution at the next time stage. 
/// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
/// @param[in, out] kkt_residual Split KKT reisdual at the current time stage. 
///
template <typename ConfigVectorType, typename TangentVectorType>
void linearizeBackwardEuler(
    const double dt, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const SplitSolution& s, const SplitSolution& s_next, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual);

///
/// @brief Linearizes the state equation of backward Euler at the terminal stage. 
/// @param[in] dt Time step. 
/// @param[in] q_prev Configuration at the previous time stage. 
/// @param[in] v_prev Generalized velocity at the previous time stage. 
/// @param[in] s Solution at the current tiem stage. 
/// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
/// @param[in, out] kkt_residual Split KKT reisdual at the current time stage. 
///
template <typename ConfigVectorType, typename TangentVectorType>
void linearizeBackwardEulerTerminal(
    const double dt, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual);

///
/// @brief Computes the residual in the state equation of forward Euler. 
/// @param[in] dt Time step. 
/// @param[in] s Solution at the current time stage. 
/// @param[in] q_next Configuration at the next time stage. 
/// @param[in] v_next Generalized velocity at the next time stage. 
/// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
///
template <typename ConfigVectorType, typename TangentVectorType>
void computeForwardEulerResidual(
    const double dt, const SplitSolution& s, 
    const Eigen::MatrixBase<ConfigVectorType>& q_next, 
    const Eigen::MatrixBase<TangentVectorType>& v_next, 
    SplitKKTResidual& kkt_residual);

///
/// @brief Computes the residual in the state equation of backward Euler. 
/// @param[in] dt Time step. 
/// @param[in] q_prev Configuration at the previous time stage. 
/// @param[in] v_prev Generalized velocity at the previous time stage. 
/// @param[in] s Solution at the current time stage. 
/// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
///
template <typename ConfigVectorType, typename TangentVectorType>
void computeBackwardEulerResidual(
    const double dt, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual);

} // namespace stateequation 
} // namespace unconstr
} // namespace idocp 

#include "idocp/unconstr/unconstr_state_equation.hxx"

#endif // IDOCP_UNCONSTR_STATE_EQUATION_HPP_ 