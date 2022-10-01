#ifndef ROBOTOC_STO_COST_FUNCTION_HPP_
#define ROBOTOC_STO_COST_FUNCTION_HPP_

#include <vector>
#include <memory>

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
  /// @brief Append a cost function component to the cost function.
  /// @param[in] cost shared pointer to the switching tmie cost function  
  /// component appended to the cost.
  ///
  void push_back(const STOCostFunctionComponentBasePtr& cost);

  ///
  /// @brief Clear cost function by removing all components.
  ///
  void clear();

  ///
  /// @brief Computes the cost on the switching times. 
  /// @param[in] time_discretization Time discretization of the hybrid optimal 
  /// control problem.
  /// @return Cost on the switching times.
  ///
  double evalCost(const TimeDiscretization& time_discretization) const;

  ///
  /// @brief Computes the cost and its derivative on the switching times. 
  /// @param[in] time_discretization Time discretization of the hybrid optimal 
  /// control problem.
  /// @param[out] lt Derivative of the cost w.r.t. the switching times.
  ///
  double linearizeCost(const TimeDiscretization& time_discretization,
                       Eigen::VectorXd& lt) const;

  ///
  /// @brief Computes the cost and its derivative on the switching times. 
  /// @param[in] time_discretization Time discretization of the hybrid optimal 
  /// control problem.
  /// @param[out] lt Derivative of the cost w.r.t. the switching times.
  /// @param[out] Qtt Hessian of the cost w.r.t. the switching times.
  ///
  double quadratizeCost(const TimeDiscretization& time_discretization,
                        Eigen::VectorXd& lt, Eigen::MatrixXd& Qtt) const;

private:
  std::vector<STOCostFunctionComponentBasePtr> costs_;
};

} // namespace robotoc

#endif // ROBOTOC_STO_COST_FUNCTION_HPP_ 