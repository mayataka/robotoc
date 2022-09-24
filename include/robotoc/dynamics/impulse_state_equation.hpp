#ifndef ROBOTOC_IMPULSE_STATE_EQUATION_HPP_
#define ROBOTOC_IMPULSE_STATE_EQUATION_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/se3_jacobian_inverse.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/dynamics/state_equation_data.hpp"
#include "robotoc/dynamics/state_equation.hpp"


namespace robotoc {

///
/// @brief Computes the residual in the impulse state equation. 
/// @param[in] robot Robot model. 
/// @param[in] s Solution at the current impulse stage. 
/// @param[in] q_next Configuration at the next time stage. 
/// @param[in] v_next Generalized velocity at the next time stage. 
/// @param[in, out] kkt_residual Impulse split KKT residual at the current 
/// impulse stage. 
///
void evalImpulseStateEquation(const Robot& robot, const SplitSolution& s, 
                              const Eigen::VectorXd& q_next, 
                              const Eigen::VectorXd& v_next, 
                              SplitKKTResidual& kkt_residual);

///
/// @brief Linearizes the impulse state equation. 
/// @param[in] robot Robot model. 
/// @param[in] q_prev Configuration at the previous time stage. 
/// @param[in] s Solution at the current impulse stage. 
/// @param[in] s_next Solution at the next time stage. 
/// @param[in, out] data Data structure for the state equation.
/// @param[in, out] kkt_matrix Impulse split KKT matrix at the current impulse 
/// stage. 
/// @param[in, out] kkt_residual Impulse split KKT residual at the current 
/// impulse stage. 
///
void linearizeImpulseStateEquation(const Robot& robot,  
                                   const Eigen::VectorXd& q_prev, 
                                   const SplitSolution& s, 
                                   const SplitSolution& s_next, 
                                   StateEquationData& data, 
                                   SplitKKTMatrix& kkt_matrix, 
                                   SplitKKTResidual& kkt_residual);

///
/// @brief Corrects the linearized state equation using the Jacobian of the 
/// Lie group. 
/// @param[in] robot Robot model. 
/// @param[in] s Solution at the current impulse stage. 
/// @param[in] s_next Solution at the next time stage. 
/// @param[in, out] data Data structure for the state equation.
/// @param[in, out] kkt_matrix Impulse split KKT matrix at the current impulse 
/// stage. 
/// @param[in, out] kkt_residual Impulse split KKT residual at the current 
/// impulse stage. 
///
void correctLinearizeImpulseStateEquation(const Robot& robot, 
                                          const SplitSolution& s, 
                                          const SplitSolution& s_next, 
                                          StateEquationData& data,
                                          SplitKKTMatrix& kkt_matrix, 
                                          SplitKKTResidual& kkt_residual);

///
/// @brief Corrects the costate direction using the Jacobian of the Lie group. 
/// @param[in, out] data Data structure for the state equation.
/// @param[in, out] d Split direction. 
///
void correctCostateDirection(StateEquationData& data, SplitDirection& d);

} // namespace robotoc 

#endif // ROBOTOC_IMPULSE_STATE_EQUATION_HPP_ 