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


namespace robotoc {

///
/// @brief Linearizes the state equation at the teminal stage. 
/// @param[in] robot Robot model. 
/// @param[in] q_prev Configuration at the previous time stage. 
/// @param[in] s Solution at the current stage. 
/// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
/// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
///
void linearizeTerminalStateEquation(const Robot& robot, StateEquationData& data, 
                                    const Eigen::VectorXd& q_prev, 
                                    const SplitSolution& s, 
                                    SplitKKTMatrix& kkt_matrix, 
                                    SplitKKTResidual& kkt_residual);

///
/// @brief Corrects the linearized state equation using the Jacobian of the 
/// Lie group. 
/// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
///
void correctLinearizeTerminalStateEquation(StateEquationData& data, 
                                           SplitKKTMatrix& kkt_matrix);

///
/// @class TerminalStateEquation
/// @brief State equation at the terminal stage. Only represent kinematic relation
/// between stages.
///
class TerminalStateEquation {
public:
  ///
  /// @brief Constructs a terminal state equation.
  /// @param[in] robot Robot model. 
  ///
  TerminalStateEquation(const Robot& robot);

  ///
  /// @brief Default constructor.  
  ///
  TerminalStateEquation();

  ///
  /// @brief Destructor. 
  ///
  ~TerminalStateEquation() = default;

  ///
  /// @brief Default copy constructor. 
  ///
  TerminalStateEquation(const TerminalStateEquation&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  TerminalStateEquation& operator=(const TerminalStateEquation&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  TerminalStateEquation(TerminalStateEquation&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  TerminalStateEquation& operator=(TerminalStateEquation&&) noexcept = default;

  ///
  /// @brief Linearizes the state equation at the teminal stage. 
  /// @param[in] robot Robot model. 
  /// @param[in] q_prev Configuration at the previous time stage. 
  /// @param[in] s Solution at the current stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
  /// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
  ///
  static void linearize(const Robot& robot, StateEquationData& data, 
                        const Eigen::VectorXd& q_prev, const SplitSolution& s, 
                        SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual);

  ///
  /// @brief Corrects the linearized state equation using the Jacobian of the 
  /// Lie group. 
  /// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
  ///
  static void correctLinearize(StateEquationData& data, SplitKKTMatrix& kkt_matrix);

  ///
  /// @brief Corrects the costate direction using the Jacobian of the Lie group. 
  /// @param[in, out] d Split direction at the terminal stage. 
  ///
  static void correctCostateDirection(StateEquationData& data, SplitDirection& d);

  ///
  /// @brief Linearizes the state equation at the teminal stage. 
  /// @param[in] robot Robot model. 
  /// @param[in] q_prev Configuration at the previous time stage. 
  /// @param[in] s Solution at the current stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
  /// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
  ///
  void linearizeStateEquation(const Robot& robot, const Eigen::VectorXd& q_prev, 
                              const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
                              SplitKKTResidual& kkt_residual);

  ///
  /// @brief Corrects the linearized state equation using the Jacobian of the 
  /// Lie group. 
  /// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
  ///
  void correctLinearizedStateEquation(SplitKKTMatrix& kkt_matrix);

  ///
  /// @brief Corrects the costate direction using the Jacobian of the Lie group. 
  /// @param[in, out] d Split direction at the terminal stage. 
  ///
  void correctCostateDirection(SplitDirection& d);

private:
  StateEquationData data_;

};

} // namespace robotoc 

#endif // ROBOTOC_TERMINAL_STATE_EQUATION_HPP_ 