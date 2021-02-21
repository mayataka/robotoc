#ifndef IDOCP_TIME_VARYING_TASK_SPACE_3D_COST_HPP_
#define IDOCP_TIME_VARYING_TASK_SPACE_3D_COST_HPP_

#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

class TimeVaryingTaskSpace3DRefBase {
public:
  TimeVaryingTaskSpace3DRefBase() {}

  virtual ~TimeVaryingTaskSpace3DRefBase() {}

  TimeVaryingTaskSpace3DRefBase(const TimeVaryingTaskSpace3DRefBase&) = default;

  TimeVaryingTaskSpace3DRefBase& operator=(
      const TimeVaryingTaskSpace3DRefBase&) = default;

  TimeVaryingTaskSpace3DRefBase(
      TimeVaryingTaskSpace3DRefBase &&) noexcept = default;

  TimeVaryingTaskSpace3DRefBase& operator=(
      TimeVaryingTaskSpace3DRefBase&&) noexcept = default;

  virtual void compute_q_3d_ref(const double t, 
                                Eigen::VectorXd& q_3d_ref) const = 0;
};


class TimeVaryingTaskSpace3DCost final : public CostFunctionComponentBase {
public:
  TimeVaryingTaskSpace3DCost(
      const Robot& robot, const int frame_id, 
      const std::shared_ptr<TimeVaryingTaskSpace3DRefBase>& ref);

  TimeVaryingTaskSpace3DCost();

  ~TimeVaryingTaskSpace3DCost();

  // Use defalut copy constructor.
  TimeVaryingTaskSpace3DCost(const TimeVaryingTaskSpace3DCost&) = default;

  // Use defalut copy operator.
  TimeVaryingTaskSpace3DCost& operator=(
      const TimeVaryingTaskSpace3DCost&) = default;

  // Use defalut move constructor.
  TimeVaryingTaskSpace3DCost(TimeVaryingTaskSpace3DCost&&) noexcept = default;

  // Use defalut copy operator.
  TimeVaryingTaskSpace3DCost& operator=(
      TimeVaryingTaskSpace3DCost&&) noexcept = default;

  bool useKinematics() const override;

  void set_ref(const std::shared_ptr<TimeVaryingTaskSpace3DRefBase>& ref);

  void set_q_3d_weight(const Eigen::Vector3d& q_3d_weight);

  void set_qf_3d_weight(const Eigen::Vector3d& qf_3d_weight);

  void set_qi_3d_weight(const Eigen::Vector3d& qi_3d_weight);

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
  std::shared_ptr<TimeVaryingTaskSpace3DRefBase> ref_;
  Eigen::Vector3d q_3d_weight_, qf_3d_weight_, qi_3d_weight_;

};

} // namespace idocp


#endif // IDOCP_TIME_VARYING_TASK_SPACE_3D_COST_HPP_ 