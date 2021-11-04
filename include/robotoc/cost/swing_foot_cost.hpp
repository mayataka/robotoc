#ifndef ROBOTOC_SWING_FOOT_COST_HPP_
#define ROBOTOC_SWING_FOOT_COST_HPP_

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


namespace robotoc {

///
/// @class SwingFootCost
/// @brief Cost on the task space position. 
///
class SwingFootCost final : public CostFunctionComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] frame_id Frame of interest.
  ///
  SwingFootCost(const Robot& robot, const int frame_id);

  ///
  /// @brief Default constructor. 
  ///
  SwingFootCost();

  ///
  /// @brief Destructor. 
  ///
  ~SwingFootCost();

  ///
  /// @brief Default copy constructor. 
  ///
  SwingFootCost(const SwingFootCost&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SwingFootCost& operator=(const SwingFootCost&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SwingFootCost(SwingFootCost&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SwingFootCost& operator=(SwingFootCost&&) noexcept = default;

  ///
  /// @brief Sets the reference position. 
  /// @param[in] q_3d_ref Reference position.
  ///
  void set_q_3d_ref(const Eigen::Vector3d& q_3d_ref);

  ///
  /// @brief Sets the weight vector. 
  /// @param[in] q_3d_weight Weight vector on the position error. 
  ///
  void set_q_weight(const Eigen::Vector3d& q_3d_weight);

  ///
  /// @brief Sets the terminal weight vector. 
  /// @param[in] qf_3d_weight Terminal weight vector on the position error. 
  ///
  void set_qf_weight(const Eigen::Vector3d& qf_3d_weight);

  ///
  /// @brief Sets the weight vector at impulse. 
  /// @param[in] qi_3d_weight Weight vector on the position error at impulse. 
  ///
  void set_qi_weight(const Eigen::Vector3d& qi_3d_weight);

  double evalStageCost(Robot& robot, const ContactStatus& contact_status, 
                       CostFunctionData& data, const double t, const double dt, 
                       const SplitSolution& s) const override;

  void evalStageCostDerivatives(Robot& robot, const ContactStatus& contact_status, 
                                CostFunctionData& data, const double t, 
                                const double dt, const SplitSolution& s, 
                                SplitKKTResidual& kkt_residual) const override;

  void evalStageCostHessian(Robot& robot, const ContactStatus& contact_status, 
                            CostFunctionData& data, const double t, 
                            const double dt, const SplitSolution& s, 
                            SplitKKTMatrix& kkt_matrix) const override;

  double evalTerminalCost(Robot& robot, CostFunctionData& data, 
                          const double t, const SplitSolution& s) const override;

  void evalTerminalCostDerivatives(Robot& robot, CostFunctionData& data, 
                                   const double t, const SplitSolution& s, 
                                   SplitKKTResidual& kkt_residual) const override;

  void evalTerminalCostHessian(Robot& robot, CostFunctionData& data, 
                               const double t, const SplitSolution& s, 
                               SplitKKTMatrix& kkt_matrix) const override;

  double evalImpulseCost(Robot& robot, const ImpulseStatus& impulse_status, 
                         CostFunctionData& data, const double t, 
                         const ImpulseSplitSolution& s) const override;

  void evalImpulseCostDerivatives(Robot& robot, const ImpulseStatus& impulse_status, 
                                  CostFunctionData& data, const double t, 
                                  const ImpulseSplitSolution& s, 
                                  ImpulseSplitKKTResidual& kkt_residual) const;

  void evalImpulseCostHessian(Robot& robot, const ImpulseStatus& impulse_status, 
                              CostFunctionData& data, const double t, 
                              const ImpulseSplitSolution& s, 
                              ImpulseSplitKKTMatrix& kkt_matrix) const override;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int frame_id_;
  Eigen::Vector3d q_3d_ref_, q_3d_weight_, qf_3d_weight_, qi_3d_weight_;

};

} // namespace robotoc

#endif // ROBOTOC_SWING_FOOT_COST_HPP_