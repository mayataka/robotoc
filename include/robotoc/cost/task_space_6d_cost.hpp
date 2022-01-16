#ifndef ROBOTOC_TASK_SPACE_6D_COST_HPP_
#define ROBOTOC_TASK_SPACE_6D_COST_HPP_

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
#include "robotoc/impulse/impulse_split_kkt_residual.hpp"
#include "robotoc/impulse/impulse_split_kkt_matrix.hpp"


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
  /// @param[in] trans_ref Reference translation.
  /// @param[in] rot_ref Reference rotation matrix.
  ///
  void set_x6d_ref(const Eigen::Vector3d& trans_ref, 
                   const Eigen::Matrix3d& rot_ref);

  ///
  /// @brief Sets the weight vectors. 
  /// @param[in] trans_weight Weight vector on the position error. 
  /// @param[in] rot_weight Weight vector on the rotation error. 
  ///
  void set_x6d_weight(const Eigen::Vector3d& trans_weight, 
                      const Eigen::Vector3d& rot_weight);

  ///
  /// @brief Sets the weight vectors at the terminal stage. 
  /// @param[in] trans_weight Weight vector on the position error at the 
  /// terminal stage. 
  /// @param[in] rot_weight Weight vector on the rotation error at the 
  /// terminal stage.
  ///
  void set_x6df_weight(const Eigen::Vector3d& trans_weight, 
                       const Eigen::Vector3d& rot_weight);

  ///
  /// @brief Sets the weight vectors at the impulse stages. 
  /// @param[in] trans_weight Weight vector on the position error at the 
  /// impulse stages. 
  /// @param[in] rot_weight Weight vector on the rotation error at the 
  /// impulse stages.
  ///
  void set_x6di_weight(const Eigen::Vector3d& trans_weight, 
                       const Eigen::Vector3d& rot_weight);

  bool useKinematics() const override;

  double evalStageCost(Robot& robot, const ContactStatus& contact_status, 
                       CostFunctionData& data, const int time_stage_in_phase, 
                       const double t, const double dt, 
                       const SplitSolution& s) const override;

  void evalStageCostDerivatives(Robot& robot, const ContactStatus& contact_status, 
                                CostFunctionData& data, const int time_stage_in_phase, 
                                const double t, const double dt, 
                                const SplitSolution& s, 
                                SplitKKTResidual& kkt_residual) const override;

  void evalStageCostHessian(Robot& robot, const ContactStatus& contact_status, 
                            CostFunctionData& data, const int time_stage_in_phase, 
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
  SE3 x6d_ref_, x6d_ref_inv_;
  Eigen::VectorXd x6d_weight_, x6df_weight_, x6di_weight_;

};

} // namespace robotoc


#endif // ROBOTOC_TASK_SPACE_6D_COST_HPP_