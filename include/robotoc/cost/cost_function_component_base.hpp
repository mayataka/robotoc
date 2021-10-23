#ifndef ROBOTOC_COST_FUNCTION_COMPONENT_BASE_HPP_
#define ROBOTOC_COST_FUNCTION_COMPONENT_BASE_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/impulse/impulse_split_solution.hpp"
#include "robotoc/impulse/impulse_split_kkt_residual.hpp"
#include "robotoc/impulse/impulse_split_kkt_matrix.hpp"


namespace robotoc {

///
/// @class CostFunctionComponentBase
/// @brief Base class of components of cost function.
///
class CostFunctionComponentBase {
public:
  ///
  /// @brief Default constructor. 
  ///
  CostFunctionComponentBase() {}

  ///
  /// @brief Destructor. 
  ///
  virtual ~CostFunctionComponentBase() {}

  ///
  /// @brief Default copy constructor. 
  ///
  CostFunctionComponentBase(const CostFunctionComponentBase&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  CostFunctionComponentBase& operator=(const CostFunctionComponentBase&) 
      = default;

  ///
  /// @brief Default move constructor. 
  ///
  CostFunctionComponentBase(CostFunctionComponentBase&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  CostFunctionComponentBase& operator=(CostFunctionComponentBase&&) noexcept 
      = default;

  ///
  /// @brief Check if the cost function component requres kinematics
  /// (forward kinematics and its Jacobians) of robot model.
  /// @return true if the cost function component requres kinematics of 
  /// Robot model. false if not.
  ///
  virtual bool useKinematics() const = 0;

  ///
  /// @brief Computes the stage cost. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dt Time step.
  /// @param[in] s Split solution.
  /// @return Stage cost.
  ///
  virtual double computeStageCost(Robot& robot, CostFunctionData& data, 
                                  const double t, const double dt, 
                                  const SplitSolution& s) const = 0;

  ///
  /// @brief Computes the first-order partial derivatives of the stage cost. 
  /// This function is always called just after computeStageCost().
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dt Time step.
  /// @param[in] s Split solution.
  /// @param[in, out] kkt_residual Split KKT residual. The partial derivatives 
  /// are added to this object.
  ///
  virtual void computeStageCostDerivatives(
      Robot& robot, CostFunctionData& data, const double t, const double dt, 
      const SplitSolution& s, SplitKKTResidual& kkt_residual) const = 0;

  ///
  /// @brief Computes the Hessian, i.e., the second-order partial derivatives of 
  /// the stage cost.  This function is always called just after 
  /// computeStageCostDerivatives().
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dt Time step.
  /// @param[in] s Split solution.
  /// @param[in, out] kkt_matrix Split KKT matrix. The Hessians are added to 
  /// this object.
  ///
  virtual void computeStageCostHessian(
      Robot& robot, CostFunctionData& data, const double t, const double dt, 
      const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const = 0;

  ///
  /// @brief Computes the terminal cost. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @return Terminal cost.
  ///
  virtual double computeTerminalCost(Robot& robot, CostFunctionData& data, 
                                     const double t, 
                                     const SplitSolution& s) const = 0;

  ///
  /// @brief Computes the first-order partial derivatives of the terminal cost. 
  /// This function is always called just after computeTerminalCost().
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[in, out] kkt_residual Split KKT residual. The partial derivatives 
  /// are added to this object.
  ///
  virtual void computeTerminalCostDerivatives(
      Robot& robot, CostFunctionData& data, const double t, 
      const SplitSolution& s, SplitKKTResidual& kkt_residual) const = 0;

  ///
  /// @brief Computes the Hessian, i.e., the second-order partial derivatives of 
  /// the teminal cost. This function is always called just after 
  /// computeTerminalCostDerivatives().
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[in, out] kkt_matrix Split KKT matrix. The Hessians are added to 
  /// this object.
  ///
  virtual void computeTerminalCostHessian(
      Robot& robot, CostFunctionData& data, const double t, 
      const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const = 0;

  ///
  /// @brief Computes the impulse cost. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @return Impulse cost.
  ///
  virtual double computeImpulseCost(Robot& robot, CostFunctionData& data, 
                                    const double t, 
                                    const ImpulseSplitSolution& s) const = 0;

  ///
  /// @brief Computes the first-order partial derivatives of the impulse cost. 
  /// This function is always called just after computeImpulseCost().
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[in, out] kkt_residual Split KKT residual. The partial derivatives 
  /// are added to this object.
  ///
  virtual void computeImpulseCostDerivatives(
      Robot& robot, CostFunctionData& data, const double t, 
      const ImpulseSplitSolution& s, 
      ImpulseSplitKKTResidual& kkt_residual) const = 0;

  ///
  /// @brief Computes the Hessian, i.e., the second-order partial derivatives of 
  /// the impulse cost. This function is always called just after 
  /// computeImpulseCostDerivatives().
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[in, out] kkt_matrix Impulse split KKT matrix. The Hessians are  
  /// added to this object.
  ///
  virtual void computeImpulseCostHessian(
      Robot& robot, CostFunctionData& data, const double t, 
      const ImpulseSplitSolution& s, 
      ImpulseSplitKKTMatrix& kkt_matrix) const = 0; 

};

} // namespace robotoc

#endif // ROBOTOC_COST_FUNCTION_COMPONENT_BASE_HPP_