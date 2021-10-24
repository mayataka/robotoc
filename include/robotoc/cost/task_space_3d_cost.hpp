#ifndef ROBOTOC_TASK_SPACE_3D_COST_HPP_
#define ROBOTOC_TASK_SPACE_3D_COST_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
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
/// @class TaskSpace3DCost
/// @brief Cost on the task space position. 
///
class TaskSpace3DCost final : public CostFunctionComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] frame_id Frame of interest.
  ///
  TaskSpace3DCost(const Robot& robot, const int frame_id);

  ///
  /// @brief Default constructor. 
  ///
  TaskSpace3DCost();

  ///
  /// @brief Destructor. 
  ///
  ~TaskSpace3DCost();

  ///
  /// @brief Default copy constructor. 
  ///
  TaskSpace3DCost(const TaskSpace3DCost&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  TaskSpace3DCost& operator=(const TaskSpace3DCost&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  TaskSpace3DCost(TaskSpace3DCost&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  TaskSpace3DCost& operator=(TaskSpace3DCost&&) noexcept = default;

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

  bool useKinematics() const override;

  double evalStageCost(Robot& robot, CostFunctionData& data, const double t, 
                       const double dt, const SplitSolution& s) const override;

  void evalStageCostDerivatives(Robot& robot, CostFunctionData& data, 
                                const double t, const double dt, 
                                const SplitSolution& s, 
                                SplitKKTResidual& kkt_residual) const override;

  void evalStageCostHessian(Robot& robot, CostFunctionData& data, 
                            const double t, const double dt, 
                            const SplitSolution& s, 
                            SplitKKTMatrix& kkt_matrix) const override;

  double evalTerminalCost(Robot& robot, CostFunctionData& data, 
                          const double t, const SplitSolution& s) const override;

  void evalTerminalCostDerivatives(Robot& robot, CostFunctionData& data, 
                                   const double t, const SplitSolution& s, 
                                   SplitKKTResidual& kkt_residual) const override;

  void evalTerminalCostHessian(Robot& robot, CostFunctionData& data, 
                               const double t, const SplitSolution& s, 
                               SplitKKTMatrix& kkt_matrix) const override;

  double evalImpulseCost(Robot& robot, CostFunctionData& data, 
                         const double t, 
                         const ImpulseSplitSolution& s) const override;

  void evalImpulseCostDerivatives(Robot& robot, CostFunctionData& data, 
                                  const double t, const ImpulseSplitSolution& s, 
                                  ImpulseSplitKKTResidual& kkt_residual) const;

  void evalImpulseCostHessian(Robot& robot, CostFunctionData& data, 
                              const double t, const ImpulseSplitSolution& s, 
                              ImpulseSplitKKTMatrix& kkt_matrix) const override;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int frame_id_;
  Eigen::Vector3d q_3d_ref_, q_3d_weight_, qf_3d_weight_, qi_3d_weight_;

};

} // namespace robotoc


#endif // ROBOTOC_TASK_SPACE_3D_COST_HPP_