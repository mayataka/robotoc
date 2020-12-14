#ifndef IDOCP_JOINT_SPACE_IMPULSE_COST_HPP_
#define IDOCP_JOINT_SPACE_IMPULSE_COST_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/impulse_cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"


namespace idocp {

class JointSpaceImpulseCost final : public ImpulseCostFunctionComponentBase {
public:
  JointSpaceImpulseCost(const Robot& robot);

  JointSpaceImpulseCost();

  ~JointSpaceImpulseCost();

  // Use defalut copy constructor.
  JointSpaceImpulseCost(const JointSpaceImpulseCost&) = default;

  // Use defalut copy operator.
  JointSpaceImpulseCost& operator=(const JointSpaceImpulseCost&) = default;

  // Use defalut move constructor.
  JointSpaceImpulseCost(JointSpaceImpulseCost&&) noexcept = default;

  // Use defalut move assign operator.
  JointSpaceImpulseCost& operator=(JointSpaceImpulseCost&&) noexcept = default;

  void set_q_ref(const Eigen::VectorXd& q_ref);

  void set_v_ref(const Eigen::VectorXd& v_ref);

  void set_dv_ref(const Eigen::VectorXd& dv_ref);

  void set_q_weight(const Eigen::VectorXd& q_weight);

  void set_v_weight(const Eigen::VectorXd& v_weight);

  void set_dv_weight(const Eigen::VectorXd& dv_weight);

  double l(Robot& robot, CostFunctionData& data, const double t, 
           const ImpulseSplitSolution& s) const override;

  void lq(Robot& robot, CostFunctionData& data, const double t, 
          const ImpulseSplitSolution& s, 
          ImpulseSplitKKTResidual& kkt_residual) const override;

  void lv(Robot& robot, CostFunctionData& data, const double t, 
          const ImpulseSplitSolution& s, 
          ImpulseSplitKKTResidual& kkt_residual) const override;

  void ldv(Robot& robot, CostFunctionData& data, const double t, 
           const ImpulseSplitSolution& s, 
           ImpulseSplitKKTResidual& kkt_residual) const override;

  void lf(Robot& robot, CostFunctionData& data, const double t, 
          const ImpulseSplitSolution& s, 
          ImpulseSplitKKTResidual& kkt_residual) const override {}

  void lqq(Robot& robot, CostFunctionData& data, const double t, 
           const ImpulseSplitSolution& s, 
           ImpulseSplitKKTMatrix& kkt_matrix) const override;

  void lvv(Robot& robot, CostFunctionData& data, const double t, 
           const ImpulseSplitSolution& s, 
           ImpulseSplitKKTMatrix& kkt_matrix) const override;

  void ldvdv(Robot& robot, CostFunctionData& data, const double t, 
             const ImpulseSplitSolution& s, 
             ImpulseSplitKKTMatrix& kkt_matrix) const override;

  void lff(Robot& robot, CostFunctionData& data, const double t, 
           const ImpulseSplitSolution& s, 
           ImpulseSplitKKTMatrix& kkt_matrix) const override {}

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int dimq_, dimv_;
  Eigen::VectorXd q_ref_, v_ref_, dv_ref_, q_weight_, v_weight_, dv_weight_;
};

} // namespace idocp


#endif // IDOCP_JOINT_SPACE_IMPULSE_COST_HPP_ 