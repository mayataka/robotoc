#ifndef IDOCP_IMPULSE_STATE_EQUATION_HPP_
#define IDOCP_IMPULSE_STATE_EQUATION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/lie_derivative_inverter.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"


namespace idocp {

///
/// @class ImpulseStateEquation
/// @brief State equation. Only represent kinematic relation
/// between stages.
///
class ImpulseStateEquation {
public:
  ///
  /// @brief Constructs an impulse state equation.
  /// @param[in] robot Robot model. 
  ///
  ImpulseStateEquation(const Robot& robot);

  ///
  /// @brief Default constructor.  
  ///
  ImpulseStateEquation();
  
  ///
  /// @brief Destructor. 
  ///
  ~ImpulseStateEquation();

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseStateEquation(const ImpulseStateEquation&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  ImpulseStateEquation& operator=(const ImpulseStateEquation&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseStateEquation(ImpulseStateEquation&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseStateEquation& operator=(ImpulseStateEquation&&) noexcept = default;

  ///
  /// @brief Computes the residual in the impulse state equation. 
  /// @param[in] robot Robot model. 
  /// @param[in] s Solution at the current impulse stage. 
  /// @param[in] q_next Configuration at the next time stage. 
  /// @param[in] v_next Generalized velocity at the next time stage. 
  /// @param[in, out] kkt_residual Impulse split KKT residual at the current 
  /// impulse stage. 
  ///
  template <typename ConfigVectorType, typename TangentVectorType>
  static void computeStateEquationResidual(
      const Robot& robot, const ImpulseSplitSolution& s, 
      const Eigen::MatrixBase<ConfigVectorType>& q_next, 
      const Eigen::MatrixBase<TangentVectorType>& v_next, 
      ImpulseSplitKKTResidual& kkt_residual);

  ///
  /// @brief Linearizes the impulse state equation. 
  /// @param[in] robot Robot model. 
  /// @param[in] q_prev Configuration at the previous time stage. 
  /// @param[in] s Solution at the current impulse stage. 
  /// @param[in] s_next Solution at the next time stage. 
  /// @param[in, out] kkt_matrix Impulse split KKT matrix at the current impulse 
  /// stage. 
  /// @param[in, out] kkt_residual Impulse split KKT residual at the current 
  /// impulse stage. 
  ///
  template <typename ConfigVectorType>
  static void linearizeStateEquation(
      const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
      const ImpulseSplitSolution& s, const SplitSolution& s_next, 
      ImpulseSplitKKTMatrix& kkt_matrix, 
      ImpulseSplitKKTResidual& kkt_residual);

  ///
  /// @brief Linearizes the impulse state equation, and multiplies the inverse
  /// of the Lie derivative to ImpulseSplitKKTMatrix::Fqq and 
  /// ImpulseSplitKKTResidual::Fq.
  /// @param[in] robot Robot model. 
  /// @param[in] q_prev Configuration at the previous time stage. 
  /// @param[in] s Solution at the current impulse stage. 
  /// @param[in] s_next Solution at the next time stage. 
  /// @param[in, out] kkt_matrix Impulse split KKT matrix at the current impulse 
  /// stage. 
  /// @param[in, out] kkt_residual Impulse split KKT residual at the current 
  /// impulse stage. 
  ///
  template <typename ConfigVectorType>
  void linearizeStateEquationAlongLieGroup(
      const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
      const ImpulseSplitSolution& s, const SplitSolution& s_next, 
      ImpulseSplitKKTMatrix& kkt_matrix, 
      ImpulseSplitKKTResidual& kkt_residual);

  ///
  /// @brief Corrects the costate direction using the derivatives of the Lie group. 
  /// @param[in, out] d Split direction. 
  ///
  void correctCostateDirection(ImpulseSplitDirection& d);

private:
  Eigen::MatrixXd Fqq_inv_, Fqq_prev_inv_, Fqq_tmp_;  
  Eigen::VectorXd Fq_tmp_;
  LieDerivativeInverter lie_der_inverter_;
  bool has_floating_base_;

};

} // namespace idocp 

#include "idocp/impulse/impulse_state_equation.hxx"

#endif // IDOCP_IMPULSE_STATE_EQUATION_HPP_ 