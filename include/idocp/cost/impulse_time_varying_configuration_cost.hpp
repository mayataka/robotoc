#ifndef IDOCP_IMPULSE_TIME_VARYING_CONFIGURATION_COST_HPP_
#define IDOCP_IMPULSE_TIME_VARYING_CONFIGURATION_COST_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/impulse_cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"


namespace idocp {

class ImpulseTimeVaryingConfigurationCost final : public ImpulseCostFunctionComponentBase {
public:
  ImpulseTimeVaryingConfigurationCost(const Robot& robot);

  ImpulseTimeVaryingConfigurationCost();

  ~ImpulseTimeVaryingConfigurationCost();

  // Use defalut copy constructor.
  ImpulseTimeVaryingConfigurationCost(
      const ImpulseTimeVaryingConfigurationCost&) = default;

  // Use defalut copy operator.
  ImpulseTimeVaryingConfigurationCost& operator=(
      const ImpulseTimeVaryingConfigurationCost&) = default;

  // Use defalut move constructor.
  ImpulseTimeVaryingConfigurationCost(
      ImpulseTimeVaryingConfigurationCost&&) noexcept = default;

  // Use defalut move assign operator.
  ImpulseTimeVaryingConfigurationCost& operator=(
      ImpulseTimeVaryingConfigurationCost&&) noexcept = default;

  void set_ref(const double t0, const Eigen::VectorXd q0, 
               const Eigen::VectorXd v0);

  void set_q_weight(const Eigen::VectorXd& q_weight);

  void set_v_weight(const Eigen::VectorXd& v_weight);

  void set_dv_weight(const Eigen::VectorXd& dv_weight);

  double l(Robot& robot, CostFunctionData& data, const double t, 
           const ImpulseSplitSolution& s) const override;

  void lq(Robot& robot, CostFunctionData& data, const double t, 
          const ImpulseSplitSolution& s, 
          ImpulseKKTResidual& kkt_residual) const override;

  void lv(Robot& robot, CostFunctionData& data, const double t, 
          const ImpulseSplitSolution& s, 
          ImpulseKKTResidual& kkt_residual) const override;

  void ldv(Robot& robot, CostFunctionData& data, const double t, 
           const ImpulseSplitSolution& s, 
           ImpulseKKTResidual& kkt_residual) const override;

  void lf(Robot& robot, CostFunctionData& data, const double t, 
          const ImpulseSplitSolution& s, 
          ImpulseKKTResidual& kkt_residual) const override {}

  void lqq(Robot& robot, CostFunctionData& data, const double t, 
           const ImpulseSplitSolution& s, 
           ImpulseKKTMatrix& kkt_matrix) const override;

  void lvv(Robot& robot, CostFunctionData& data, const double t, 
           const ImpulseSplitSolution& s, 
           ImpulseKKTMatrix& kkt_matrix) const override;

  void ldvdv(Robot& robot, CostFunctionData& data, const double t, 
             const ImpulseSplitSolution& s, 
             ImpulseKKTMatrix& kkt_matrix) const override;

  void lff(Robot& robot, CostFunctionData& data, const double t, 
           const ImpulseSplitSolution& s, 
           ImpulseKKTMatrix& kkt_matrix) const override {}

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int dimq_, dimv_;
  double t0_;
  Eigen::VectorXd q0_, v0_, q_weight_, v_weight_, dv_weight_;

};

} // namespace idocp

#endif // IDOCP_IMPULSE_TIME_VARYING_CONFIGURATION_COST_HPP_ 