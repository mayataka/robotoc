#ifndef IDOCP_COST_FUNCTION_HPP_
#define IDOCP_COST_FUNCTION_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"


namespace idocp {

///
/// @class CostFunction
/// @brief Stack of the cost function. Composed by cost function components 
/// that inherits CostFunctionComponentBase.
///
class CostFunction {
public:
  using CostFunctionComponentBasePtr 
      = std::shared_ptr<CostFunctionComponentBase>;

  ///
  /// @brief Default constructor. 
  ///
  CostFunction();

  ///
  /// @brief Destructor. 
  ///
  ~CostFunction();

  ///
  /// @brief Default copy constructor. 
  ///
  CostFunction(const CostFunction&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  CostFunction& operator=(const CostFunction&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  CostFunction(CostFunction&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  CostFunction& operator=(CostFunction&&) noexcept = default;

  ///
  /// @brief Append a cost function component to the cost function.
  /// @param[in] cost shared pointer to the cost function component appended 
  /// to the cost.
  ///
  void push_back(const CostFunctionComponentBasePtr& cost);

  ///
  /// @brief Clear cost function by removing all components.
  ///
  void clear();

  ///
  /// @brief Check if the cost function component requres kinematics of robot 
  /// model.
  /// @return true if the cost function component requres kinematics of 
  /// Robot model. false if not.
  ///
  bool useKinematics() const;

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
  /// @param[in] dtau Time step.
  /// @param[in] s Split solution.
  /// @return Stage cost.
  ///
  double computeStageCost(Robot& robot, CostFunctionData& data, const double t, 
                          const double dtau, const SplitSolution& s) const;

  ///
  /// @brief Computes and returns terminal cost. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @return Terminal cost.
  ///
  double computeTerminalCost(Robot& robot, CostFunctionData& data,
                             const double t, const SplitSolution& s) const;

  ///
  /// @brief Computes and returns stage cost. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @return Stage cost.
  ///
  double computeImpulseCost(Robot& robot, CostFunctionData& data, 
                            const double t, 
                            const ImpulseSplitSolution& s) const;

  ///
  /// @brief Computes the partial derivatives of the stage cost with respect
  /// to the configuration, velocity, acceleration, and contact forces. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dtau Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkt_residual The KKT residual. The partial derivatives are 
  /// added to this data.
  ///
  void computeStageCostDerivatives(Robot& robot, CostFunctionData& data, 
                                   const double t, const double dtau, 
                                   const SplitSolution& s, 
                                   SplitKKTResidual& kkt_residual) const;

  ///
  /// @brief Computes the partial derivatives of the terminal cost with respect
  /// to the configuration and velocity. 
  /// @param[in] robot Robot modol.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[out] kkt_residual The KKT residual. The partial derivatives are 
  /// added to this data.
  ///
  void computeTerminalCostDerivatives(Robot& robot, CostFunctionData& data, 
                                      const double t, const SplitSolution& s, 
                                      SplitKKTResidual& kkt_residual) const;


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
  void computeImpulseCostDerivatives(Robot& robot, CostFunctionData& data, 
                                     const double t, 
                                     const ImpulseSplitSolution& s, 
                                     ImpulseSplitKKTResidual& kkt_residual) const;


  ///
  /// @brief Computes the Hessians of the stage cost with respect
  /// to the configuration, velocity, acceleration, and contact forces. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dtau Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkt_matrix The KKT matrix. The Hessians are added to this 
  /// data.
  ///
  void computeStageCostHessian(Robot& robot, CostFunctionData& data, 
                               const double t, const double dtau, 
                               const SplitSolution& s, 
                               SplitKKTMatrix& kkt_matrix) const;

  ///
  /// @brief Computes the Hessians of the terminal cost with respect
  /// to the configuration and velocity. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[out] kkt_matrix The KKT matrix. The Hessians are added to this 
  /// data.
  ///
  void computeTerminalCostHessian(Robot& robot, CostFunctionData& data, 
                                  const double t, const SplitSolution& s, 
                                  SplitKKTMatrix& kkt_matrix) const;

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
                               ImpulseSplitKKTMatrix& kkt_matrix) const;

  ///
  /// @brief Computes the Hessians of the stage cost with respect
  /// to the configuration. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[out] kkt_matrix The KKT matrix. The Hessians are added to this 
  /// data.
  ///
  void computeImpulseCostHessian(Robot& robot, CostFunctionData& data, 
                                 const double t, const ImpulseSplitSolution& s, 
                                 ImpulseSplitKKTMatrix& kkt_matrix) const; 

private:
  std::vector<CostFunctionComponentBasePtr> costs_;

};

} // namespace idocp

#include "idocp/cost/cost_function.hxx"

#endif // IDOCP_COST_FUNCTION_HPP_