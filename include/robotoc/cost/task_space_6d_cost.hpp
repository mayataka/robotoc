#ifndef ROBOTOC_TASK_SPACE_6D_COST_HPP_
#define ROBOTOC_TASK_SPACE_6D_COST_HPP_

#include <string>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/robot/se3.hpp"
#include "robotoc/cost/cost_function_component_base.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/impulse/impulse_split_solution.hpp"
#include "robotoc/cost/task_space_6d_ref_base.hpp"


namespace robotoc {

///
/// @class TaskSpace6DCost
/// @brief Cost on the task space placement (SE(3)). 
///
class TaskSpace6DCost final : public CostFunctionComponentBase {
public:
  using Vector6d = Eigen::Matrix<double, 6, 1>;

  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] frame_id Frame of interest.
  ///
  TaskSpace6DCost(const Robot& robot, const int frame_id);

  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] frame_name Frame of interest.
  ///
  TaskSpace6DCost(const Robot& robot, const std::string& frame_name);
///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] frame_id Frame of interest.
  /// @param[in] ref Reference task-space placement.
  ///
  TaskSpace6DCost(const Robot& robot, const int frame_id,
                  const std::shared_ptr<TaskSpace6DRefBase>& ref);

  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] frame_id Frame of interest.
  /// @param[in] const_ref Const reference task-space placement.
  ///
  TaskSpace6DCost(const Robot& robot, const int frame_id,
                  const SE3& const_ref);

  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] frame_id Frame of interest.
  /// @param[in] const_position_ref Const reference task-space position.
  /// @param[in] const_rotation_ref Const reference task-space rotation.
  ///
  TaskSpace6DCost(const Robot& robot, const int frame_id,
                  const Eigen::Vector3d& const_position_ref,
                  const Eigen::Matrix3d& const_rotation_ref);

  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] frame_name Frame of interest.
  /// @param[in] ref Reference task-space placement.
  ///
  TaskSpace6DCost(const Robot& robot, const std::string& frame_name,
                  const std::shared_ptr<TaskSpace6DRefBase>& ref);

  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] frame_name Frame of interest.
  /// @param[in] const_ref Const reference task-space placement.
  ///
  TaskSpace6DCost(const Robot& robot, const std::string& frame_name,
                  const SE3& const_ref);

  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] frame_name Frame of interest.
  /// @param[in] const_position_ref Const reference task-space position.
  /// @param[in] const_rotation_ref Const reference task-space rotation.
  ///
  TaskSpace6DCost(const Robot& robot, const std::string& frame_name,
                  const Eigen::Vector3d& const_position_ref,
                  const Eigen::Matrix3d& const_rotation_ref);

  ///
  /// @brief Default constructor. 
  ///
  TaskSpace6DCost();

  ///
  /// @brief Destructor. 
  ///
  ~TaskSpace6DCost();

  ///
  /// @brief Default copy constructor. 
  ///
  TaskSpace6DCost(const TaskSpace6DCost&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  TaskSpace6DCost& operator=(const TaskSpace6DCost&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  TaskSpace6DCost(TaskSpace6DCost&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  TaskSpace6DCost& operator=(TaskSpace6DCost&&) noexcept = default;

  ///
  /// @brief Sets the reference task-space placement. 
  /// @param[in] ref Reference task-space placement.
  ///
  void set_ref(const std::shared_ptr<TaskSpace6DRefBase>& ref);

  ///
  /// @brief Sets the const reference task-space placement. 
  /// @param[in] const_ref Const reference task-space placement.
  ///
  void set_const_ref(const SE3& const_ref);

  ///
  /// @brief Sets the const reference task-space placement. 
  /// @param[in] const_position_ref Const reference task-space position.
  /// @param[in] const_rotation_ref Const reference task-space rotation.
  ///
  void set_const_ref(const Eigen::Vector3d& const_position_ref,
                     const Eigen::Matrix3d& const_rotation_ref);

  ///
  /// @brief Sets the weight vector. 
  /// @param[in] weight_position Weight vector on the task-space position error. 
  /// @param[in] weight_rotation Weight vector on the task-space rotation error. 
  ///
  void set_weight(const Eigen::Vector3d& weight_position,
                  const Eigen::Vector3d& weight_rotation);

  ///
  /// @brief Sets the weight vector at the terminal stage. 
  /// @param[in] weight_position_terminal Weight vector on the task-space 
  /// position error at the terminal stage. 
  /// @param[in] weight_rotation_terminal Weight vector on the task-space 
  /// rotation error at the terminal stage. 
  ///
  void set_weight_terminal(const Eigen::Vector3d& weight_position_terminal,
                           const Eigen::Vector3d& weight_rotation_terminal);

  ///
  /// @brief Sets the weight vector at the impulse stage. 
  /// @param[in] weight_position_impulse Weight vector on the task-space 
  /// position error at the impulse stage. 
  /// @param[in] weight_rotation_impulse Weight vector on the task-space 
  /// rotation error at the impulse stage. 
  ///
  void set_weight_impulse(const Eigen::Vector3d& weight_position_impulse,
                          const Eigen::Vector3d& weight_rotation_impulse);

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
        ref_->updateRef(grid_info, data.x6d_ref);
        data.x6d_ref_inv = data.x6d_ref.inverse();
        data.diff_x6d = data.x6d_ref_inv * robot.framePlacement(frame_id_);
        data.diff_6d = Log6Map(data.diff_x6d);
      }
    }
    else {
      data.diff_x6d = const_ref_inv_ * robot.framePlacement(frame_id_);
      data.diff_6d = Log6Map(data.diff_x6d);
    }
  }

  double evalStageCost(Robot& robot, const ContactStatus& contact_status, 
                       CostFunctionData& data, const GridInfo& grid_info, 
                       const SplitSolution& s) const override;

  void evalStageCostDerivatives(Robot& robot, const ContactStatus& contact_status, 
                                CostFunctionData& data, const GridInfo& grid_info,
                                const SplitSolution& s, 
                                SplitKKTResidual& kkt_residual) const override;

  void evalStageCostHessian(Robot& robot, const ContactStatus& contact_status, 
                            CostFunctionData& data, const GridInfo& grid_info,  
                            const SplitSolution& s, 
                            SplitKKTMatrix& kkt_matrix) const override;

  double evalTerminalCost(Robot& robot, CostFunctionData& data, 
                          const GridInfo& grid_info, 
                          const SplitSolution& s) const override;

  void evalTerminalCostDerivatives(Robot& robot, CostFunctionData& data, 
                                   const GridInfo& grid_info, 
                                   const SplitSolution& s, 
                                   SplitKKTResidual& kkt_residual) const override;

  void evalTerminalCostHessian(Robot& robot, CostFunctionData& data, 
                               const GridInfo& grid_info, 
                               const SplitSolution& s, 
                               SplitKKTMatrix& kkt_matrix) const override;

  double evalImpulseCost(Robot& robot, const ImpulseStatus& impulse_status, 
                         CostFunctionData& data, const GridInfo& grid_info, 
                         const ImpulseSplitSolution& s) const override;

  void evalImpulseCostDerivatives(Robot& robot, const ImpulseStatus& impulse_status, 
                                  CostFunctionData& data, const GridInfo& grid_info,
                                  const ImpulseSplitSolution& s, 
                                  SplitKKTResidual& kkt_residual) const override;

  void evalImpulseCostHessian(Robot& robot, const ImpulseStatus& impulse_status, 
                              CostFunctionData& data, const GridInfo& grid_info,
                              const ImpulseSplitSolution& s, 
                              SplitKKTMatrix& kkt_matrix) const override;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int frame_id_;
  SE3 const_ref_, const_ref_inv_;
  Eigen::VectorXd weight_, weight_terminal_, weight_impulse_;
  std::shared_ptr<TaskSpace6DRefBase> ref_;
  bool use_nonconst_ref_, enable_cost_, enable_cost_terminal_, enable_cost_impulse_;
};

} // namespace robotoc


#endif // ROBOTOC_TASK_SPACE_6D_COST_HPP_