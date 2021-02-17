#ifndef IDOCP_CONFIGURATION_SPACE_COST_HPP_
#define IDOCP_CONFIGURATION_SPACE_COST_HPP_

#include "Eigen/Core"

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

class ConfigurationSpaceCost final : public CostFunctionComponentBase {
public:
  using Vector6d = Eigen::Matrix<double, 6, 1>;

  ConfigurationSpaceCost(const Robot& robot);

  ConfigurationSpaceCost();

  ~ConfigurationSpaceCost();

  // Use defalut copy constructor.
  ConfigurationSpaceCost(const ConfigurationSpaceCost&) = default;

  // Use defalut copy operator.
  ConfigurationSpaceCost& operator=(const ConfigurationSpaceCost&) = default;

  // Use defalut move constructor.
  ConfigurationSpaceCost(ConfigurationSpaceCost&&) noexcept = default;

  // Use defalut move assign operator.
  ConfigurationSpaceCost& operator=(ConfigurationSpaceCost&&) noexcept = default;

  bool useKinematics() const override;

  void set_q_ref(const Eigen::VectorXd& q_ref);

  void set_v_ref(const Eigen::VectorXd& v_ref);

  void set_u_ref(const Eigen::VectorXd& u_ref);

  void set_q_weight(const Eigen::VectorXd& q_weight);

  void set_v_weight(const Eigen::VectorXd& v_weight);

  void set_a_weight(const Eigen::VectorXd& a_weight);

  void set_u_weight(const Eigen::VectorXd& u_weight);

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

private:
  int dimq_, dimv_, dimu_;
  Eigen::VectorXd q_ref_, v_ref_, u_ref_, q_weight_, v_weight_, a_weight_, 
                  u_weight_, qf_weight_, vf_weight_, qi_weight_, vi_weight_, 
                  dvi_weight_;
};

} // namespace idocp


#endif // IDOCP_CONFIGURATION_SPACE_COST_HPP_ 