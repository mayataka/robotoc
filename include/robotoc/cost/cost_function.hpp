#ifndef ROBOTOC_COST_FUNCTION_HPP_
#define ROBOTOC_COST_FUNCTION_HPP_

#include <vector>
#include <memory>
#include <cmath>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/cost/cost_function_component_base.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/impulse/impulse_split_solution.hpp"
#include "robotoc/impulse/impulse_split_kkt_residual.hpp"
#include "robotoc/impulse/impulse_split_kkt_matrix.hpp"
#include "robotoc/hybrid/grid_info.hpp"


namespace robotoc {

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
  /// @brief Constructor with discount factor. 
  /// @param[in] discount_factor Discount factor. Must be positive and smaller 
  /// than 1.0.
  /// @param[in] discount_time_step The cost is reduced by discount_factor as 
  /// the time proceeds to this value. Must be positive.
  ///
  CostFunction(const double discount_factor, const double discount_time_step);

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
  /// @brief Clones this to a shared ptr. 
  ///
  std::shared_ptr<CostFunction> clone() const;

  ///
  /// @brief Sets the discount facor. 
  /// @param[in] discount_factor Discount factor. Must be positive and smaller 
  /// than 1.0. Otherwise, the discounted cost is disabled.
  /// @param[in] discount_time_step The cost is reduced by discount_factor as 
  /// the time proceeds to this value. Must be positive. Otherwise, the discounted
  /// cost is disabled.
  ///
  void setDiscountFactor(const double discount_factor, 
                         const double discount_time_step);

  ///
  /// @brief Gets the discount facor. 
  /// @return The discount facor. 
  ///
  double discountFactor() const;

  ///
  /// @brief Gets the discount time step. 
  /// @return The discount time step. 
  ///
  double discountTimeStep() const;

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
  /// @brief Check if the cost function component requires kinematics 
  /// (forward kinematics and its Jacobians) of robot model.
  /// @return true if the cost function component requires kinematics of 
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
  /// @brief Computes the stage cost. 
  /// @param[in] robot Robot model.
  /// @param[in] contact_status Contact status.
  /// @param[in] data Cost function data.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @return Stage cost.
  ///
  double evalStageCost(Robot& robot, const ContactStatus& contact_status,
                       CostFunctionData& data, const GridInfo& grid_info, 
                       const SplitSolution& s) const;

  ///
  /// @brief Computes the stage cost and its first-order partial derivatives. 
  /// @param[in] robot Robot model.
  /// @param[in] contact_status Contact status.
  /// @param[in] data Cost function data.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @param[in, out] kkt_residual Split KKT residual. The partial derivatives 
  /// are added to this object.
  /// @return Stage cost.
  ///
  double linearizeStageCost(Robot& robot, const ContactStatus& contact_status, 
                            CostFunctionData& data, const GridInfo& grid_info, 
                            const SplitSolution& s, 
                            SplitKKTResidual& kkt_residual) const;

  ///
  /// @brief Computes the stage cost, its first-order partial derivatives, and 
  /// its Hessian, i.e., its second-order partial derivatives. 
  /// @param[in] robot Robot model.
  /// @param[in] contact_status Contact status.
  /// @param[in] data Cost function data.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @param[in, out] kkt_residual Split KKT residual. The partial derivatives 
  /// are added to this object.
  /// @param[in, out] kkt_matrix Split KKT matrix. The Hessians are added to 
  /// this object.
  /// @return Stage cost.
  ///
  double quadratizeStageCost(Robot& robot, const ContactStatus& contact_status, 
                             CostFunctionData& data, const GridInfo& grid_info, 
                             const SplitSolution& s, 
                             SplitKKTResidual& kkt_residual,
                             SplitKKTMatrix& kkt_matrix) const;

  ///
  /// @brief Computes the terminal cost. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @return Terminal cost.
  ///
  double evalTerminalCost(Robot& robot, CostFunctionData& data, 
                          const GridInfo& grid_info, 
                          const SplitSolution& s) const;

  ///
  /// @brief Computes the terminal cost and its first-order partial derivatives. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @param[in, out] kkt_residual Split KKT residual. The partial derivatives 
  /// are added to this object.
  /// @return Stage cost.
  ///
  double linearizeTerminalCost(Robot& robot, CostFunctionData& data, 
                               const GridInfo& grid_info, const SplitSolution& s, 
                               SplitKKTResidual& kkt_residual) const;

  ///
  /// @brief Computes the terminal cost, its first-order partial derivatives, 
  /// and its Hessian, i.e., its second-order partial derivatives. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @param[in, out] kkt_residual Split KKT residual. The partial derivatives 
  /// are added to this object.
  /// @param[in, out] kkt_matrix Split KKT matrix. The Hessians are added to 
  /// this object.
  /// @return Stage cost.
  ///
  double quadratizeTerminalCost(Robot& robot, CostFunctionData& data, 
                                const GridInfo& grid_info, const SplitSolution& s, 
                                SplitKKTResidual& kkt_residual,
                                SplitKKTMatrix& kkt_matrix) const;

  ///
  /// @brief Computes the impulse cost. 
  /// @param[in] robot Robot model.
  /// @param[in] impulse_status Impulse status.
  /// @param[in] data Cost function data.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @return Stage cost.
  ///
  double evalImpulseCost(Robot& robot, const ImpulseStatus& impulse_status, 
                         CostFunctionData& data, const GridInfo& grid_info,
                         const ImpulseSplitSolution& s) const;

  ///
  /// @brief Computes the impulse cost and its first-order partial derivatives. 
  /// @param[in] robot Robot model.
  /// @param[in] impulse_status Impulse status.
  /// @param[in] data Cost function data.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @param[in, out] kkt_residual Split KKT residual. The partial derivatives 
  /// are added to this object.
  /// @return Stage cost.
  ///
  double linearizeImpulseCost(Robot& robot, const ImpulseStatus& impulse_status, 
                              CostFunctionData& data, const GridInfo& grid_info,
                              const ImpulseSplitSolution& s, 
                              ImpulseSplitKKTResidual& kkt_residual) const;

  ///
  /// @brief Computes the impulse cost, its first-order partial derivatives, 
  /// and its Hessian, i.e., its second-order partial derivatives. 
  /// @param[in] robot Robot model.
  /// @param[in] impulse_status Impulse status.
  /// @param[in] data Cost function data.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @param[in, out] kkt_residual Split KKT residual. The partial derivatives 
  /// are added to this object.
  /// @param[in, out] kkt_matrix Split KKT matrix. The Hessians are added to 
  /// this object.
  /// @return Stage cost.
  ///
  double quadratizeImpulseCost(Robot& robot, const ImpulseStatus& impulse_status, 
                               CostFunctionData& data, const GridInfo& grid_info,
                               const ImpulseSplitSolution& s, 
                               ImpulseSplitKKTResidual& kkt_residual,
                               ImpulseSplitKKTMatrix& kkt_matrix) const;

private:
  std::vector<CostFunctionComponentBasePtr> costs_;

  double discount(const double t0, const double t) const {
    assert(t >= t0);
    return std::pow(discount_factor_, ((t-t0)/discount_time_step_));
  }

  double discount_factor_, discount_time_step_;
  bool discounted_cost_;
};

} // namespace robotoc

#endif // ROBOTOC_COST_FUNCTION_HPP_