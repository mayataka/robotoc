#ifndef IDOCP_TASK_SPACE_6D_COST_HPP_
#define IDOCP_TASK_SPACE_6D_COST_HPP_

#include "Eigen/Core"
#include "pinocchio/spatial/se3.hpp"

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
/// @class TaskSpace6DCost
/// @brief Quadratic cost on the task space pose (SE(3). 
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
  /// @brief Sets the reference pose. 
  /// @param[in] position_ref Reference position.
  /// @param[in] rotation_mat_ref Reference rotation matrix.
  ///
  void set_q_6d_ref(const Eigen::Vector3d& position_ref, 
                    const Eigen::Matrix3d& rotation_mat_ref);

  ///
  /// @brief Sets the weight vectors. 
  /// @param[in] position_weight Weight vector on the position error. 
  /// @param[in] rotation_weight Weight vector on the rotation error. 
  ///
  void set_q_6d_weight(const Eigen::Vector3d& position_weight, 
                       const Eigen::Vector3d& rotation_weight);

  ///
  /// @brief Sets the terminal weight vectors. 
  /// @param[in] position_weight Temrinal weight vector on the position error. 
  /// @param[in] rotation_weight Temrinal weight vector on the rotation error. 
  ///
  void set_qf_6d_weight(const Eigen::Vector3d& position_weight, 
                        const Eigen::Vector3d& rotation_weight);

  ///
  /// @brief Sets the weight vectors at impulse. 
  /// @param[in] position_weight Weight vector on the position error at impulse. 
  /// @param[in] rotation_weight Weight vector on the rotation error at impulse. 
  ///
  void set_qi_6d_weight(const Eigen::Vector3d& position_weight, 
                        const Eigen::Vector3d& rotation_weight);

  bool useKinematics() const override;

  double computeStageCost(Robot& robot, CostFunctionData& data, const double t, 
                          const double dt, 
                          const SplitSolution& s) const override;

  double computeTerminalCost(Robot& robot, CostFunctionData& data, 
                             const double t, 
                             const SplitSolution& s) const override;

  double computeImpulseCost(Robot& robot, CostFunctionData& data, 
                            const double t, 
                            const ImpulseSplitSolution& s) const override;

  void computeStageCostDerivatives(
      Robot& robot, CostFunctionData& data, const double t, const double dt, 
      const SplitSolution& s, SplitKKTResidual& kkt_residual) const override;

  void computeTerminalCostDerivatives(
      Robot& robot, CostFunctionData& data, const double t, 
      const SplitSolution& s, SplitKKTResidual& kkt_residual) const override;

  void computeImpulseCostDerivatives(
      Robot& robot, CostFunctionData& data, const double t, 
      const ImpulseSplitSolution& s, 
      ImpulseSplitKKTResidual& kkt_residual) const;

  void computeStageCostHessian(Robot& robot, CostFunctionData& data, 
                               const double t, const double dt, 
                               const SplitSolution& s, 
                               SplitKKTMatrix& kkt_matrix) const override;

  void computeTerminalCostHessian(Robot& robot, CostFunctionData& data, 
                                  const double t, const SplitSolution& s, 
                                  SplitKKTMatrix& kkt_matrix) const override;

  void computeImpulseCostHessian(
      Robot& robot, CostFunctionData& data, const double t, 
      const ImpulseSplitSolution& s, 
      ImpulseSplitKKTMatrix& kkt_matrix) const override;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int frame_id_;
  pinocchio::SE3 SE3_ref_, SE3_ref_inv_;
  Eigen::VectorXd q_6d_weight_, qf_6d_weight_, qi_6d_weight_;

};

} // namespace idocp


#endif // IDOCP_TASK_SPACE_6D_COST_HPP_