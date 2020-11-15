#ifndef IDOCP_IMPULSE_COST_FUNCTION_HPP_
#define IDOCP_IMPULSE_COST_FUNCTION_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/impulse_cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"


namespace idocp {

///
/// @class ImpulseCostFunction
/// @brief Stack of the cost function. Composed by cost function components 
/// that inherits ImpulseCostFunctionComponentBase.
///
class ImpulseCostFunction {
public:

  ///
  /// @brief Default constructor. 
  ///
  ImpulseCostFunction();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseCostFunction();

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseCostFunction(const ImpulseCostFunction&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseCostFunction& operator=(const ImpulseCostFunction&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseCostFunction(ImpulseCostFunction&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseCostFunction& operator=(ImpulseCostFunction&&) noexcept = default;

  ///
  /// @brief Append a cost function component to the cost function.
  /// @param[in] cost shared pointer to the cost function component appended 
  /// to the cost.
  ///
  void push_back(const std::shared_ptr<ImpulseCostFunctionComponentBase>& cost);

  ///
  /// @brief Clear cost function by removing all components.
  ///
  void clear();

  ///
  /// @brief Check whether the cost function is empty or not.
  /// @return true if the cost function is empty. false if not.
  ///
  bool isEmpty() const;

  ///
  /// @brief Creates CostFunctionData according to robot model and cost 
  /// function components. 
  /// @param[in] robot robot model.
  /// @return Cost function data.
  ///
  CostFunctionData createCostFunctionData(const Robot& robot) const;

  ///
  /// @brief Computes and returns stage cost. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @return Stage cost.
  ///
  double l(Robot& robot, CostFunctionData& data, const double t, 
           const ImpulseSplitSolution& s) const;

  ///
  /// @brief Computes the partial derivatives of the stage cost with respect
  /// to the configuration, velocity, acceleration, and contact forces. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[out] kkt_residual The KKT residual. The partial derivatives are 
  /// added to this data.
  ///
  void computeStageCostDerivatives(Robot& robot, CostFunctionData& data, 
                                   const double t, 
                                   const ImpulseSplitSolution& s, 
                                   ImpulseKKTResidual& kkt_residual) const;

  ///
  /// @brief Computes the Hessians of the stage cost with respect
  /// to the configuration, velocity, acceleration, and contact forces. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[out] kkt_matrix The KKT matrix. The Hessians are added to this 
  /// data.
  ///
  void computeStageCostHessian(Robot& robot, CostFunctionData& data, 
                               const double t, const ImpulseSplitSolution& s, 
                               ImpulseKKTMatrix& kkt_matrix) const;

private:
  std::vector<std::shared_ptr<ImpulseCostFunctionComponentBase>> costs_;

};

} // namespace idocp

#include "idocp/cost/impulse_cost_function.hxx"

#endif // IDOCP_IMPULSE_COST_FUNCTION_HPP_ 