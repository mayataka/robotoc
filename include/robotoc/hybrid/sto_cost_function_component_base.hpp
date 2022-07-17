#ifndef ROBOTOC_STO_COST_FUNCTION_COMPONENT_BASE_HPP_
#define ROBOTOC_STO_COST_FUNCTION_COMPONENT_BASE_HPP_

#include "robotoc/hybrid/time_discretization.hpp"

#include "Eigen/Core"


namespace robotoc {

#define DEFINE_DEFAULT_CLONE_STO_COST_FUNCTION_COMPONENT(Derived) \
  std::shared_ptr<STOCostFunctionComponentBase> clone() const override { return std::make_shared<Derived>(*this); } 

///
/// @class STOCostFunctionComponentBase
/// @brief Base class of components of the cost function of the switching time
/// optimization (STO) problem.
///
class STOCostFunctionComponentBase {
public:
  ///
  /// @brief Default constructor. 
  ///
  STOCostFunctionComponentBase() {}

  ///
  /// @brief Destructor. 
  ///
  virtual ~STOCostFunctionComponentBase() {}

  ///
  /// @brief Default copy constructor. 
  ///
  STOCostFunctionComponentBase(const STOCostFunctionComponentBase&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  STOCostFunctionComponentBase& operator=(
      const STOCostFunctionComponentBase&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  STOCostFunctionComponentBase(
      STOCostFunctionComponentBase&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  STOCostFunctionComponentBase& operator=(
      STOCostFunctionComponentBase&&) noexcept = default;

  ///
  /// @brief Clones this to a shared ptr. 
  ///
  virtual std::shared_ptr<STOCostFunctionComponentBase> clone() const = 0;

  ///
  /// @brief Computes the cost on the switching times. 
  /// @param[in] discretization Time discretization of the hybrid optimal 
  /// control problem.
  /// @return Cost on the switching times.
  ///
  virtual double evalCost(const TimeDiscretization& discretization) const = 0;

  ///
  /// @brief Computes the derivative of the cost on the switching times. 
  /// @param[in] discretization Time discretization of the hybrid optimal 
  /// control problem.
  /// @param[out] lts Derivative of the cost w.r.t. the switching times.
  ///
  virtual void evalCostDerivatives(const TimeDiscretization& discretization,
                                   Eigen::VectorXd& lts) const = 0;

  ///
  /// @brief Computes the twice-time derivative (Hessian) of the cost on the 
  /// switching times. 
  /// @param[in] discretization Time discretization of the hybrid optimal 
  /// control problem.
  /// @param[out] Qts Hessian of the cost w.r.t. the switching times.
  ///
  virtual void evalCostHessian(const TimeDiscretization& discretization,
                               Eigen::MatrixXd& Qts) const = 0;

};

} // namespace robotoc

#endif // ROBOTOC_STO_COST_FUNCTION_COMPONENT_BASE_HPP_ 