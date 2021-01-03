#ifndef IDOCP_TIME_VARYING_CONFIGURATION_SPACE_COST_HPP_
#define IDOCP_TIME_VARYING_CONFIGURATION_SPACE_COST_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

class TimeVaryingConfigurationSpaceCost final : public CostFunctionComponentBase {
public:
  TimeVaryingConfigurationSpaceCost(const Robot& robot);

  TimeVaryingConfigurationSpaceCost();

  ~TimeVaryingConfigurationSpaceCost();

  // Use defalut copy constructor.
  TimeVaryingConfigurationSpaceCost(
      const TimeVaryingConfigurationSpaceCost&) = default;

  // Use defalut copy operator.
  TimeVaryingConfigurationSpaceCost& operator=(
      const TimeVaryingConfigurationSpaceCost&) = default;

  // Use defalut move constructor.
  TimeVaryingConfigurationSpaceCost(
      TimeVaryingConfigurationSpaceCost&&) noexcept = default;

  // Use defalut move assign operator.
  TimeVaryingConfigurationSpaceCost& operator=(
      TimeVaryingConfigurationSpaceCost&&) noexcept = default;

  bool useKinematics() const override;

  void set_ref(const Robot& robot, const double t_begin, const double t_end, 
               const Eigen::VectorXd q_begin, const Eigen::VectorXd v);

  void set_q_weight(const Eigen::VectorXd& q_weight);

  void set_v_weight(const Eigen::VectorXd& v_weight);

  void set_a_weight(const Eigen::VectorXd& a_weight);

  void set_qf_weight(const Eigen::VectorXd& qf_weight);

  void set_vf_weight(const Eigen::VectorXd& vf_weight);

  void set_qi_weight(const Eigen::VectorXd& qi_weight);

  void set_vi_weight(const Eigen::VectorXd& vi_weight);

  void set_dvi_weight(const Eigen::VectorXd& dvi_weight);

  double computeStageCost(Robot& robot, CostFunctionData& data, const double t, 
                          const double dtau, const SplitSolution& s) const;

  double computeTerminalCost(Robot& robot, CostFunctionData& data, 
                             const double t, const SplitSolution& s) const;

  double computeImpulseCost(Robot& robot, CostFunctionData& data, 
                            const double t, 
                            const ImpulseSplitSolution& s) const;

  void computeStageCostDerivatives(Robot& robot, CostFunctionData& data, 
                                   const double t, const double dtau, 
                                   const SplitSolution& s, 
                                   SplitKKTResidual& kkt_residual) const;

  void computeTerminalCostDerivatives(Robot& robot, CostFunctionData& data, 
                                      const double t, const SplitSolution& s, 
                                      SplitKKTResidual& kkt_residual) const;

  void computeImpulseCostDerivatives(Robot& robot, CostFunctionData& data, 
                                     const double t, 
                                     const ImpulseSplitSolution& s, 
                                     ImpulseSplitKKTResidual& kkt_residual) const;

  void computeStageCostHessian(Robot& robot, CostFunctionData& data, 
                               const double t, const double dtau, 
                               const SplitSolution& s, 
                               SplitKKTMatrix& kkt_matrix) const;

  void computeTerminalCostHessian(Robot& robot, CostFunctionData& data, 
                                  const double t, const SplitSolution& s, 
                                  SplitKKTMatrix& kkt_matrix) const;

  void computeImpulseCostHessian(Robot& robot, CostFunctionData& data, 
                                 const double t, const ImpulseSplitSolution& s, 
                                 ImpulseSplitKKTMatrix& kkt_matrix) const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW


  void set_q_ref(const Robot& robot, const double t, 
                 Eigen::VectorXd& q_ref) const {
    if (t > t_begin_ && t < t_end_) {
      robot.integrateConfiguration(q_begin_, v_ref_, t-t_begin_, q_ref);
    }
    else if (t <= t_begin_) {
      q_ref = q_begin_;
    }
    else {
      q_ref = q_end_;
    }
  }

  Eigen::VectorXd v_ref(const double t) const {
    if (t > t_begin_ && t < t_end_) {
      return v_ref_;
    }
    else {
      return Eigen::VectorXd::Zero(dimv_);
    }
  }

private:
  int dimq_, dimv_;
  double t_begin_, t_end_;
  Eigen::VectorXd q_begin_, q_end_, v_ref_, q_weight_, v_weight_, a_weight_,  
                  qf_weight_, vf_weight_, qi_weight_, vi_weight_, dvi_weight_;

};

} // namespace idocp

#endif // IDOCP_TIME_VARYING_CONFIGURATION_SPACE_COST_HPP_ 