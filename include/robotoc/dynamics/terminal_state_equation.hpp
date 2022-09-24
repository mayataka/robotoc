#ifndef ROBOTOC_TERMINAL_STATE_EQUATION_HPP_
#define ROBOTOC_TERMINAL_STATE_EQUATION_HPP_

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
/// @brief Linearizes the state equation at the teminal stage. 
/// @param[in] robot Robot model. 
/// @param[in] q_prev Configuration at the previous time stage. 
/// @param[in] s Solution at the current stage. 
/// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
/// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
///
void linearizeTerminalStateEquation(const Robot& robot,  
                                    const Eigen::VectorXd& q_prev, 
                                    const SplitSolution& s, 
                                    StateEquationData& data,
                                    SplitKKTMatrix& kkt_matrix, 
                                    SplitKKTResidual& kkt_residual);

///
/// @brief Corrects the linearized state equation using the Jacobian of the 
/// Lie group. 
/// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
///
void correctLinearizeTerminalStateEquation(StateEquationData& data, 
                                           SplitKKTMatrix& kkt_matrix);

} // namespace robotoc 

#endif // ROBOTOC_TERMINAL_STATE_EQUATION_HPP_ 