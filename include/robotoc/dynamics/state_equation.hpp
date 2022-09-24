#ifndef ROBOTOC_STATE_EQUATION_HPP_
#define ROBOTOC_STATE_EQUATION_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/se3_jacobian_inverse.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/dynamics/state_equation_data.hpp"


namespace robotoc {

///
/// @brief Computes the residual in the state equation. 
/// @param[in] robot Robot model. 
/// @param[in] dt Time step. 
/// @param[in] s Solution at the current time stage. 
/// @param[in] q_next Configuration at the next time stage. 
/// @param[in] v_next Generalized velocity at the next time stage. 
/// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
///
void evalStateEquation(const Robot& robot, const double dt, 
                       const SplitSolution& s, 
                       const Eigen::VectorXd& q_next, 
                       const Eigen::VectorXd& v_next, 
                       SplitKKTResidual& kkt_residual);

///
/// @brief Linearizes the state equation. 
/// @param[in] robot Robot model. 
/// @param[in] dt Time step. 
/// @param[in] q_prev Configuration at the previous time stage. 
/// @param[in] s Solution at the current stage. 
/// @param[in] s_next Solution at the next time stage. 
/// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
/// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
///
void linearizeStateEquation(const Robot& robot, const double dt, 
                            const Eigen::VectorXd& q_prev, 
                            const SplitSolution& s, const SplitSolution& s_next, 
                            StateEquationData& data, SplitKKTMatrix& kkt_matrix, 
                            SplitKKTResidual& kkt_residual);

///
/// @brief Corrects the linearized state equation using the Jacobian of the 
/// Lie group. 
/// @param[in] robot Robot model. 
/// @param[in] dt Time step. 
/// @param[in] s Solution at the current stage. 
/// @param[in] s_next Solution at the next time stage. 
/// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
/// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
///
void correctLinearizeStateEquation(const Robot& robot, const double dt, 
                                   const SplitSolution& s, 
                                   const SplitSolution& s_next, 
                                   StateEquationData& data, 
                                   SplitKKTMatrix& kkt_matrix, 
                                   SplitKKTResidual& kkt_residual);

///
/// @brief Corrects the costate direction using the Jacobian of the Lie group. 
/// @param[in, out] d Split direction. 
///
void correctCostateDirection(StateEquationData& data, SplitDirection& d);

///
/// @brief Computes the initial state direction using the result of  
/// linearizeStateEquation() and correctLinearizeStateEquation().
/// @param[in] robot Robot model. 
/// @param[in] q0 Initial configuration. 
/// @param[in] v0 Initial generalized velocity. 
/// @param[in] s0 Split solution at the initial stage. 
/// @param[in, out] d0 Split direction at the initial stage. 
///
void computeInitialStateDirection(const Robot& robot, 
                                  const Eigen::VectorXd& q0, 
                                  const Eigen::VectorXd& v0, 
                                  const SplitSolution& s0, 
                                  const StateEquationData& data,
                                  SplitDirection& d0);

} // namespace robotoc 

#endif // ROBOTOC_STATE_EQUATION_HPP_