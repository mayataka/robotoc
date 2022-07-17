#ifndef ROBOTOC_STO_COST_FUNCTION_HPP_
#define ROBOTOC_STO_COST_FUNCTION_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "robotoc/hybrid/time_discretization.hpp"
#include "robotoc/hybrid/sto_cost_function_component_base.hpp"
#include "robotoc/ocp/kkt_residual.hpp"
#include "robotoc/ocp/kkt_matrix.hpp"


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
  /// @brief Destructor. 
  ///
  ~STOCostFunction();

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
  /// @brief Clones this to a shared ptr. 
  ///
  std::shared_ptr<STOCostFunction> clone() const;

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
  /// @param[in] discretization Time discretization of the hybrid optimal 
  /// control problem.
  /// @return Cost on the switching times.
  ///
  double evalCost(const TimeDiscretization& discretization);

  ///
  /// @brief Computes the cost on the switching times and its first-order 
  /// partial derivatives. 
  /// @param[in] discretization Time discretization of the hybrid optimal 
  /// control problem.
  /// @param[in, out] kkt_residual KKT residual. The partial derivatives 
  /// are added to this object.
  /// @return Cost on the switching times.
  ///
  double linearizeCost(const TimeDiscretization& discretization,
                       KKTResidual& kkt_residual); 

  ///
  /// @brief Computes the cost, its first-order partial derivatives, and 
  /// its Hessian, i.e., its second-order partial derivatives. 
  /// @param[in] discretization Time discretization of the hybrid optimal 
  /// control problem.
  /// @param[out] kkt_matrix KKT matrix.
  /// @param[out] kkt_residual KKT residual.
  /// @return Cost on the switching times.
  ///
  double quadratizeCost(const TimeDiscretization& discretization,
                        KKTMatrix& kkt_matrix, KKTResidual& kkt_residual);

private:
  std::vector<STOCostFunctionComponentBasePtr> costs_;
  Eigen::VectorXd lts_;
  Eigen::MatrixXd Qts_;

  void setToKKT(const TimeDiscretization& discretization,
                KKTResidual& kkt_residual);

  void setToKKT(const TimeDiscretization& discretization,
                KKTMatrix& kkt_matrix, KKTResidual& kkt_residual);

};

} // namespace robotoc

#endif // ROBOTOC_STO_COST_FUNCTION_HPP_ 