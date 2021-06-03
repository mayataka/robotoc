#ifndef IDOCP_TERMINAL_STATE_EQUATION_HPP_
#define IDOCP_TERMINAL_STATE_EQUATION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/lie_derivative_inverter.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

///
/// @class TerminalStateEquation
/// @brief State equation of forward Euler. Only represent kinematic relation
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
  /// @brief Linearizes the state equation of forward Euler at the teminal stage. 
  /// @param[in] robot Robot model. 
  /// @param[in] q_prev Configuration at the previous time stage. 
  /// @param[in] s Solution at the current stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
  /// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
  ///
  template <typename ConfigVectorType>
  static void linearizeForwardEuler(
      const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
      const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
      SplitKKTResidual& kkt_residual);

  ///
  /// @brief Linearizes the state equation of forward Euler at the teminal stage
  /// and computes the inverse matrix of the Lie derivative.
  /// @param[in] robot Robot model. 
  /// @param[in] q_prev Configuration at the previous time stage. 
  /// @param[in] s Solution at the current stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
  /// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
  ///
  template <typename ConfigVectorType>
  void linearizeForwardEulerLieDerivative(
      const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
      const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
      SplitKKTResidual& kkt_residual);

  ///
  /// @brief Corrects the costate direction using the derivatives of the Lie group. 
  /// @param[in, out] d Split direction at the terminal stage. 
  ///
  void correctCostateDirection(SplitDirection& d);

private:
  Eigen::MatrixXd Fqq_inv_, Fqq_prev_inv_;  
  Eigen::VectorXd Fq_tmp_;
  LieDerivativeInverter lie_der_inverter_;
  bool has_floating_base_;

};

} // namespace idocp 

#include "idocp/ocp/terminal_state_equation.hxx"

#endif // IDOCP_TERMINAL_STATE_EQUATION_HPP_ 