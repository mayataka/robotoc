#ifndef ROBOTOC_TASK_SPACE_3D_COST_HPP_
#define ROBOTOC_TASK_SPACE_3D_COST_HPP_

#include <string>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/impact_status.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/cost/cost_function_component_base.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/cost/task_space_3d_ref_base.hpp"


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
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] frame_name Frame of interest.
  ///
  TaskSpace3DCost(const Robot& robot, const std::string& frame_name);

  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] frame_id Frame of interest.
  /// @param[in] ref Reference task-space position.
  ///
  TaskSpace3DCost(const Robot& robot, const int frame_id, 
                  const std::shared_ptr<TaskSpace3DRefBase>& ref);

  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] frame_id Frame of interest.
  /// @param[in] const_ref Const reference task-space position.
  ///
  TaskSpace3DCost(const Robot& robot, const int frame_id, 
                  const Eigen::Vector3d& const_ref);

  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] frame_name Frame of interest.
  /// @param[in] ref Reference task-space position.
  ///
  TaskSpace3DCost(const Robot& robot, const std::string& frame_name,
                  const std::shared_ptr<TaskSpace3DRefBase>& ref);

  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] frame_name Frame of interest.
  /// @param[in] const_ref Const reference task-space position.
  ///
  TaskSpace3DCost(const Robot& robot, const std::string& frame_name,
                  const Eigen::Vector3d& const_ref);

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
  /// @brief Sets the reference task-space position. 
  /// @param[in] ref Reference task-space position.
  ///
  void set_ref(const std::shared_ptr<TaskSpace3DRefBase>& ref);

  ///
  /// @brief Sets the const reference task-space position. 
  /// @param[in] const_ref Const reference task-space position.
  ///
  void set_const_ref(const Eigen::Vector3d& const_ref);

  ///
  /// @brief Sets the weight vector. 
  /// @param[in] weight Weight vector on the task-space position error. 
  ///
  void set_weight(const Eigen::Vector3d& weight);

  ///
  /// @brief Sets the weight vector for the terminal stage. 
  /// @param[in] weight_terminal Weight vector on the task-space position error 
  /// at the terminal stage. 
  ///
  void set_weight_terminal(const Eigen::Vector3d& weight_terminal);

  ///
  /// @brief Sets the weight vector for the impact stage. 
  /// @param[in] weight_impact Weight vector on the task-space position error 
  /// at the impact stage. 
  ///
  void set_weight_impact(const Eigen::Vector3d& weight_impact);

  ///
  /// @brief Evaluate if the cost is active for given grid_info. 
  /// @param[in] grid_info Grid info.
  /// @return Cost status (if the cost is active or not).
  ///
  bool isCostActive(const GridInfo& grid_info) const {
    if (use_nonconst_ref_) {
      return ref_->isActive(grid_info);
    }
    else {
      return true;
    }
  }

  ///
  /// @brief Evaluate the difference between the robot's task-space position 
  /// status and reference. 
  /// @param[in] robot Robot model.
  /// @param[in, out] data Cost funciton data.
  /// @param[in] grid_info Grid info
  ///
  void evalDiff(const Robot& robot, CostFunctionData& data, 
                const GridInfo& grid_info) const {
    if (use_nonconst_ref_) {
      if (ref_->isActive(grid_info)) {
        ref_->updateRef(grid_info, data.x3d_ref);
        data.diff_3d = robot.framePosition(frame_id_) - data.x3d_ref;
      }
    }
    else {
      data.diff_3d = robot.framePosition(frame_id_) - const_ref_;
    }
  }

  double evalStageCost(Robot& robot, const ContactStatus& contact_status, 
                       const GridInfo& grid_info, const SplitSolution& s,
                       CostFunctionData& data) const override;

  void evalStageCostDerivatives(Robot& robot, const ContactStatus& contact_status, 
                                const GridInfo& grid_info, const SplitSolution& s, 
                                CostFunctionData& data,
                                SplitKKTResidual& kkt_residual) const override;

  void evalStageCostHessian(Robot& robot, const ContactStatus& contact_status, 
                            const GridInfo& grid_info, const SplitSolution& s, 
                            CostFunctionData& data,
                            SplitKKTMatrix& kkt_matrix) const override;

  double evalTerminalCost(Robot& robot, const GridInfo& grid_info,
                          const SplitSolution& s,
                          CostFunctionData& data) const override;

  void evalTerminalCostDerivatives(Robot& robot, const GridInfo& grid_info,
                                   const SplitSolution& s, CostFunctionData& data, 
                                   SplitKKTResidual& kkt_residual) const override;

  void evalTerminalCostHessian(Robot& robot, const GridInfo& grid_info,
                               const SplitSolution& s, CostFunctionData& data, 
                               SplitKKTMatrix& kkt_matrix) const override;

  double evalImpactCost(Robot& robot, const ImpactStatus& impact_status, 
                         const GridInfo& grid_info, const SplitSolution& s,
                         CostFunctionData& data) const override;

  void evalImpactCostDerivatives(Robot& robot, const ImpactStatus& impact_status, 
                                 const GridInfo& grid_info, const SplitSolution& s, 
                                 CostFunctionData& data,
                                 SplitKKTResidual& kkt_residual) const override;

  void evalImpactCostHessian(Robot& robot, const ImpactStatus& impact_status, 
                             const GridInfo& grid_info, const SplitSolution& s, 
                             CostFunctionData& data,
                             SplitKKTMatrix& kkt_matrix) const override;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int frame_id_;
  Eigen::Vector3d const_ref_, weight_, weight_terminal_, weight_impact_;
  std::shared_ptr<TaskSpace3DRefBase> ref_;
  bool use_nonconst_ref_, enable_cost_, enable_cost_terminal_, enable_cost_impact_;
};

} // namespace robotoc


#endif // ROBOTOC_TASK_SPACE_3D_COST_HPP_