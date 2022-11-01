#ifndef ROBOTOC_STO_COST_FUNCTION_COMPONENT_BASE_HPP_
#define ROBOTOC_STO_COST_FUNCTION_COMPONENT_BASE_HPP_

#include <memory>
#include <stdexcept>

#include "Eigen/Core"

#include "robotoc/ocp/time_discretization.hpp"


namespace robotoc {

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
  /// @brief Computes the cost on the switching times. 
  /// @param[in] time_discretization Time discretization of the hybrid optimal 
  /// control problem.
  /// @return Cost on the switching times.
  ///
  virtual double evalCost(const TimeDiscretization& time_discretization) const = 0;

  ///
  /// @brief Computes the derivative of the cost on the switching times. 
  /// @param[in] time_discretization Time discretization of the hybrid optimal 
  /// control problem.
  /// @param[in, out] lt The derivatives of the Lagrangian w.r.t. the switching
  /// times. 
  ///
  virtual void evalCostDerivatives(const TimeDiscretization& time_discretization,
                                   Eigen::VectorXd& lt) const = 0;

  ///
  /// @brief Computes the twice-time derivative (Hessian) of the cost on the 
  /// switching times. 
  /// @param[in] time_discretization Time discretization of the hybrid optimal 
  /// control problem.
  /// @param[in, out] Qtt The Hessian of the Lagrangian w.r.t. the switching
  /// times. 
  ///
  virtual void evalCostHessian(const TimeDiscretization& time_discretization,
                               Eigen::MatrixXd& Qtt) const = 0;

  ///
  /// @brief Gets the shared ptr of this object as the specified type. If this 
  /// fails in dynamic casting, throws an exception.
  /// @tparam Derived The derived type.
  /// @return shared ptr of this object as the specified type. 
  ///
  template <typename Derived>
  std::shared_ptr<Derived> as() const {
    Derived* derived_ptr = dynamic_cast<Derived*>(this);
    if (derived_ptr == nullptr) {
        throw std::runtime_error("[STOCostFunctionComponentBase] runtime error: failed in down-casting!");
    }
    return std::shared_ptr<Derived>(derived_ptr);
  }

};

} // namespace robotoc

#endif // ROBOTOC_STO_COST_FUNCTION_COMPONENT_BASE_HPP_ 