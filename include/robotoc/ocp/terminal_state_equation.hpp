#ifndef ROBOTOC_TERMINAL_STATE_EQUATION_HPP_
#define ROBOTOC_TERMINAL_STATE_EQUATION_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/se3_jacobian_inverse.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"


namespace robotoc {

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
  ~TerminalStateEquation();

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
  template <typename ConfigVectorType>
  static void linearizeStateEquation(
      const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
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
  Eigen::MatrixXd Fqq_inv_, Fqq_prev_inv_;  
  Eigen::VectorXd Fq_tmp_;
  SE3JacobianInverse se3_jac_inverse_;
  bool has_floating_base_;

};

} // namespace robotoc 

#include "robotoc/ocp/terminal_state_equation.hxx"

#endif // ROBOTOC_TERMINAL_STATE_EQUATION_HPP_ 