#ifndef ROBOTOC_STO_COST_FUNCTION_HPP_
#define ROBOTOC_STO_COST_FUNCTION_HPP_

#include <vector>
#include <memory>
#include <unordered_map>

#include "Eigen/Core"

#include "robotoc/ocp/time_discretization.hpp"
#include "robotoc/sto/sto_cost_function_component_base.hpp"
#include "robotoc/core/kkt_residual.hpp"
#include "robotoc/core/kkt_matrix.hpp"


namespace robotoc {

///
/// @class STOCostFunction
/// @brief Stack of the cost function of the switching time optimization (STO)
///  problem. Composed by cost function components that inherits 
/// STOCostFunctionComponentBase.
///
class STOCostFunction {
public:
  using STOCostFunctionComponentBasePtr 
      = std::shared_ptr<STOCostFunctionComponentBase>;

  ///
  /// @brief Default constructor. 
  ///
  STOCostFunction();

  ///
  /// @brief Default destructor. 
  ///
  ~STOCostFunction() = default;

  ///
  /// @brief Default copy constructor. 
  ///
  STOCostFunction(const STOCostFunction&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  STOCostFunction& operator=(const STOCostFunction&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  STOCostFunction(STOCostFunction&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  STOCostFunction& operator=(STOCostFunction&&) noexcept = default;

  ///
  /// @brief Checks if thsi has a STO cost function component of the specified 
  /// name. 
  /// @param[in] name Name of the STO cost function component.
  /// @return treu if a STO cost function component of the specified name exists. 
  ///
  bool exist(const std::string& name) const;

  ///
  /// @brief Adds a STO cost function component. If a component of the same name 
  /// exists, throws an exeption.
  /// @param[in] name Name of the STO cost function component.
  /// @param[in] cost shared pointer to the STO cost function component.
  ///
  void add(const std::string& name, const STOCostFunctionComponentBasePtr& cost);

  ///
  /// @brief Erases a STO cost function component. If a component of the 
  /// specified name does not exist, throws an exeption.
  /// @param[in] name Name of the STO cost function component.
  ///
  void erase(const std::string& name);

  ///
  /// @brief Gets a STO cost function component. If a component of the specified 
  /// name does not exist, throws an exeption. 
  /// @param[in] name Name of the STO cost function component.
  /// @return Shared ptr to the specified STO cost function component.
  ///
  STOCostFunctionComponentBasePtr get(const std::string& name) const;

  ///
  /// @brief Clear cost function by removing all components.
  ///
  void clear();

  ///
  /// @brief Computes the cost on the switching times. 
  /// @param[in] time_discretization Time discretization.
  /// @return Cost on the switching times.
  ///
  double evalCost(const TimeDiscretization& time_discretization) const;

  ///
  /// @brief Computes the cost and its derivative on the switching times. 
  /// @param[in] time_discretization Time discretization.
  /// @param[in, out] lt The derivatives of the Lagrangian w.r.t. the switching
  /// times. 
  ///
  double linearizeCost(const TimeDiscretization& time_discretization,
                       Eigen::VectorXd& lt) const;

  ///
  /// @brief Computes the cost and its derivative on the switching times. 
  /// @param[in] time_discretization Time discretization.
  /// @param[in, out] lt The derivatives of the Lagrangian w.r.t. the switching
  /// times. 
  /// @param[in, out] Qtt The Hessian of the Lagrangian w.r.t. the switching
  /// times. 
  ///
  double quadratizeCost(const TimeDiscretization& time_discretization,
                        Eigen::VectorXd& lt, Eigen::MatrixXd& Qtt) const;

private:
  std::vector<STOCostFunctionComponentBasePtr> costs_;
  std::unordered_map<std::string, size_t> cost_names_;

};

} // namespace robotoc

#endif // ROBOTOC_STO_COST_FUNCTION_HPP_ 