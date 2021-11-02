#ifndef ROBOTOC_SWITCHING_TIME_COST_FUNCTION_COMPONENT_BASE_HPP_
#define ROBOTOC_SWITCHING_TIME_COST_FUNCTION_COMPONENT_BASE_HPP_

#include "robotoc/hybrid/hybrid_ocp_discretization.hpp"

#include "Eigen/Core"


namespace robotoc {

///
/// @class SwitchingTimeCostFunctionComponentBase
/// @brief Base class of components of the cost function on the switching time.
///
class SwitchingTimeCostFunctionComponentBase {
public:
  ///
  /// @brief Default constructor. 
  ///
  SwitchingTimeCostFunctionComponentBase() {}

  ///
  /// @brief Destructor. 
  ///
  virtual ~SwitchingTimeCostFunctionComponentBase() {}

  ///
  /// @brief Default copy constructor. 
  ///
  SwitchingTimeCostFunctionComponentBase(
      const SwitchingTimeCostFunctionComponentBase&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SwitchingTimeCostFunctionComponentBase& operator=(
      const SwitchingTimeCostFunctionComponentBase&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SwitchingTimeCostFunctionComponentBase(
      SwitchingTimeCostFunctionComponentBase&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SwitchingTimeCostFunctionComponentBase& operator=(
      SwitchingTimeCostFunctionComponentBase&&) noexcept = default;

  ///
  /// @brief Computes the cost on the switching times. 
  /// @param[in] discretization Discretization of the optimal control problem.
  /// @return Cost on the switching times.
  ///
  virtual double evalCost(
      const HybridOCPDiscretization& discretization) const = 0;

  ///
  /// @brief Computes the derivative of the cost on the switching times. 
  /// @param[in] discretization Discretization of the optimal control problem.
  /// @param[out] lts Derivative of the cost w.r.t. the switching times.
  ///
  virtual void evalCostDerivatives(const HybridOCPDiscretization& discretization,
                                   Eigen::VectorXd& lts) const = 0;

  ///
  /// @brief Computes the twice-time derivative (Hessian) of the cost on the 
  /// switching times. 
  /// @param[in] discretization Discretization of the optimal control problem.
  /// @param[out] Qts Hessian of the cost w.r.t. the switching times.
  ///
  virtual void evalCostHessian(const HybridOCPDiscretization& discretization,
                               Eigen::MatrixXd& Qts) const = 0;

};

} // namespace robotoc

#endif // ROBOTOC_SWITCHING_TIME_COST_FUNCTION_COMPONENT_BASE_HPP_ 