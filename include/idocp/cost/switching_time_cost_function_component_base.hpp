#ifndef IDOCP_SWITCHING_TIME_COST_FUNCTION_COMPONENT_BASE_HPP_
#define IDOCP_SWITCHING_TIME_COST_FUNCTION_COMPONENT_BASE_HPP_

#include "Eigen/Core"


namespace idocp {

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
  /// @param[in] ts Switching times.
  /// @return Cost on the switching times.
  ///
  virtual double computeCost(const double t0, const double tf, 
                             const Eigen::VectorXd& ts) const = 0;

  ///
  /// @brief Computes the derivative of the cost on the switching times. 
  /// @param[in] ts Switching times.
  /// @param[in] hts Derivative of the cost w.r.t. the switching times.
  ///
  virtual void computeCostDerivatives(const double t0, const double tf, 
                                      const Eigen::VectorXd& ts,
                                      Eigen::VectorXd& hts) const = 0;

  ///
  /// @brief Computes the twice-time derivative (Hessian) of the cost on the 
  /// switching times. 
  /// @param[in] ts Switching times.
  /// @param[in] Qts Hessian of the cost w.r.t. the switching times.
  ///
  virtual void computeCostHessian(const double t0, const double tf, 
                                  const Eigen::VectorXd& ts,
                                  Eigen::MatrixXd& Qts) const = 0;

};

} // namespace idocp

#endif // IDOCP_SWITCHING_TIME_COST_FUNCTION_COMPONENT_BASE_HPP_ 