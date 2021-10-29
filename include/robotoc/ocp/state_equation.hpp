#ifndef ROBOTOC_STATE_EQUATION_HPP_
#define ROBOTOC_STATE_EQUATION_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/lie_derivative_inverter.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"


namespace robotoc {

///
/// @class StateEquation
/// @brief State equation. Only represent kinematic relation between stages.
///
class StateEquation {
public:
  ///
  /// @brief Constructs a state equation.
  /// @param[in] robot Robot model. 
  ///
  StateEquation(const Robot& robot);

  ///
  /// @brief Default constructor.  
  ///
  StateEquation();
  
  ///
  /// @brief Destructor. 
  ///
  ~StateEquation();

  ///
  /// @brief Default copy constructor. 
  ///
  StateEquation(const StateEquation&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  StateEquation& operator=(const StateEquation&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  StateEquation(StateEquation&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  StateEquation& operator=(StateEquation&&) noexcept = default;

  ///
  /// @brief Computes the residual in the state equation. 
  /// @param[in] robot Robot model. 
  /// @param[in] dt Time step. 
  /// @param[in] s Solution at the current time stage. 
  /// @param[in] q_next Configuration at the next time stage. 
  /// @param[in] v_next Generalized velocity at the next time stage. 
  /// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
  ///
  template <typename ConfigVectorType, typename TangentVectorType>
  static void evalStateEquation(
      const Robot& robot, const double dt, const SplitSolution& s, 
      const Eigen::MatrixBase<ConfigVectorType>& q_next, 
      const Eigen::MatrixBase<TangentVectorType>& v_next, 
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
  template <typename ConfigVectorType, typename SplitSolutionType>
  static void linearizeStateEquation(
      const Robot& robot, const double dt, 
      const Eigen::MatrixBase<ConfigVectorType>& q_prev, const SplitSolution& s, 
      const SplitSolutionType& s_next, SplitKKTMatrix& kkt_matrix, 
      SplitKKTResidual& kkt_residual);

  ///
  /// @brief Linearizes the state equation, and multiplies the inverse of the
  /// Lie derivative to SplitKKTMatrix::Fqq, SplitKKTMatrix::Fqv, and 
  /// SplitKKTResidual::Fq.
  /// @param[in] robot Robot model. 
  /// @param[in] dt Time step. 
  /// @param[in] q_prev Configuration at the previous time stage. 
  /// @param[in] s Solution at the current stage. 
  /// @param[in] s_next Solution at the next time stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
  /// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
  ///
  template <typename ConfigVectorType, typename SplitSolutionType>
  void linearizeStateEquationAlongLieGroup(
      const Robot& robot, const double dt, 
      const Eigen::MatrixBase<ConfigVectorType>& q_prev, const SplitSolution& s, 
      const SplitSolutionType& s_next, SplitKKTMatrix& kkt_matrix, 
      SplitKKTResidual& kkt_residual);

  ///
  /// @brief Corrects the costate direction using the derivatives of the Lie group. 
  /// @param[in, out] d Split direction. 
  ///
  void correctCostateDirection(SplitDirection& d);

  ///
  /// @brief Computes the initial state direction using the result of  
  /// StateEquation::linearizeStateEquationAlongLieGroup().
  /// @param[in] robot Robot model. 
  /// @param[in] q0 Initial configuration. 
  /// @param[in] v0 Initial generalized velocity. 
  /// @param[in] s0 Split solution at the initial stage. 
  /// @param[in, out] d0 Split direction at the initial stage. 
  ///
  template <typename ConfigVectorType, typename TangentVectorType>
  void computeInitialStateDirection(
      const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q0, 
      const Eigen::MatrixBase<TangentVectorType>& v0, 
      const SplitSolution& s0, SplitDirection& d0) const;

private:
  Eigen::MatrixXd Fqq_inv_, Fqq_prev_inv_, Fqq_tmp_;  
  Eigen::VectorXd Fq_tmp_;
  LieDerivativeInverter lie_der_inverter_;
  bool has_floating_base_;

};

} // namespace robotoc 

#include "robotoc/ocp/state_equation.hxx"

#endif // ROBOTOC_STATE_EQUATION_HPP_