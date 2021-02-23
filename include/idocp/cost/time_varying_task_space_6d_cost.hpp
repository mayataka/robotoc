#ifndef IDOCP_TIME_VARYING_TASK_SPACE_6D_COST_HPP_
#define IDOCP_TIME_VARYING_TASK_SPACE_6D_COST_HPP_

#include <memory>

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

class TimeVaryingTaskSpace6DRefBase {
public:
  TimeVaryingTaskSpace6DRefBase() {}

  virtual ~TimeVaryingTaskSpace6DRefBase() {}

  TimeVaryingTaskSpace6DRefBase(const TimeVaryingTaskSpace6DRefBase&) = default;

  TimeVaryingTaskSpace6DRefBase& operator=(
      const TimeVaryingTaskSpace6DRefBase&) = default;

  TimeVaryingTaskSpace6DRefBase(
      TimeVaryingTaskSpace6DRefBase &&) noexcept = default;

  TimeVaryingTaskSpace6DRefBase& operator=(
      TimeVaryingTaskSpace6DRefBase&&) noexcept = default;

  virtual void compute_q_6d_ref(const double t, 
                                pinocchio::SE3& se3_ref) const = 0;
};


class TimeVaryingTaskSpace6DCost final : public CostFunctionComponentBase {
public:
  using Vector6d = Eigen::Matrix<double, 6, 1>;

  TimeVaryingTaskSpace6DCost(
      const Robot& robot, const int frame_id,
      const std::shared_ptr<TimeVaryingTaskSpace6DRefBase>& ref);

  TimeVaryingTaskSpace6DCost();

  ~TimeVaryingTaskSpace6DCost();

  // Use defalut copy constructor.
  TimeVaryingTaskSpace6DCost(const TimeVaryingTaskSpace6DCost&) = default;

  // Use defalut copy operator.
  TimeVaryingTaskSpace6DCost& operator=(
      const TimeVaryingTaskSpace6DCost&) = default;

  // Use defalut move constructor.
  TimeVaryingTaskSpace6DCost(TimeVaryingTaskSpace6DCost&&) noexcept = default;

  // Use defalut copy operator.
  TimeVaryingTaskSpace6DCost& operator=(
      TimeVaryingTaskSpace6DCost&&) noexcept = default;

  bool useKinematics() const override;

  void set_ref(const std::shared_ptr<TimeVaryingTaskSpace6DRefBase>& ref);

  void set_q_6d_weight(const Eigen::Vector3d& position_weight, 
                       const Eigen::Vector3d& rotation_weight);

  void set_qf_6d_weight(const Eigen::Vector3d& position_weight, 
                        const Eigen::Vector3d& rotation_weight);

  void set_qi_6d_weight(const Eigen::Vector3d& position_weight, 
                        const Eigen::Vector3d& rotation_weight);

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
  std::shared_ptr<TimeVaryingTaskSpace6DRefBase> ref_;
  Eigen::VectorXd q_6d_weight_, qf_6d_weight_, qi_6d_weight_;

};

} // namespace idocp


#endif // IDOCP_TIME_VARYING_TASK_SPACE_6D_COST_HPP_ 