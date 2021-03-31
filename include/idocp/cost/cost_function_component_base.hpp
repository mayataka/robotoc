#ifndef IDOCP_COST_FUNCTION_COMPONENT_BASE_HPP_
#define IDOCP_COST_FUNCTION_COMPONENT_BASE_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"


namespace idocp {

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
  /// @brief Check if the cost function component requres kinematics of robot 
  /// model.
  /// @return true if the cost function component requres kinematics of 
  /// Robot model. false if not.
  ///
  virtual bool useKinematics() const = 0;

  ///
  /// @brief Computes and returns stage cost. 
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
  /// @brief Computes and returns terminal cost. 
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
  /// @brief Computes and returns stage cost. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @return Stage cost.
  ///
  virtual double computeImpulseCost(Robot& robot, CostFunctionData& data, 
                                    const double t, 
                                    const ImpulseSplitSolution& s) const = 0;

  ///
  /// @brief Computes the partial derivatives of the stage cost with respect
  /// to the configuration. 
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
  /// @brief Computes the partial derivatives of the terminal cost with respect
  /// to the velocity.
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
  /// @brief Computes the partial derivatives of the stage cost with respect
  /// to the velocity. 
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
  /// @brief Computes the Hessians of the stage cost with respect
  /// to the configuration. 
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
  /// @brief Computes the Hessians of the terminal cost with respect
  /// to the configuration. 
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
  /// @brief Computes the Hessians of the stage cost with respect
  /// to the configuration. 
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

} // namespace idocp

#endif // IDOCP_COST_FUNCTION_COMPONENT_BASE_HPP_