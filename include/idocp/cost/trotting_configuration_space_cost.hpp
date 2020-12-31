#ifndef IDOCP_TROTTING_CONFIGURATION_SPACE_COST_HPP_
#define IDOCP_TROTTING_CONFIGURATION_SPACE_COST_HPP_

#include <cmath>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

struct TrottingSwingAngles {
  double front_swing_thigh;
  double front_swing_knee; 
  double front_stance_thigh; 
  double front_stance_knee; 
  double hip_swing_thigh;
  double hip_swing_knee;
  double hip_stance_thigh;
  double hip_stance_knee;

  TrottingSwingAngles()
    : front_swing_thigh(0), front_swing_knee(0),
      front_stance_thigh(0), front_stance_knee(0), 
      hip_swing_thigh(0), hip_swing_knee(0),
      hip_stance_thigh(0), hip_stance_knee(0) {}
};


class TrottingConfigurationSpaceCost final : public CostFunctionComponentBase {
public:

  TrottingConfigurationSpaceCost(const Robot& robot);

  TrottingConfigurationSpaceCost();

  ~TrottingConfigurationSpaceCost();

  // Use defalut copy constructor.
  TrottingConfigurationSpaceCost(const TrottingConfigurationSpaceCost&) = default;

  // Use defalut copy operator.
  TrottingConfigurationSpaceCost& operator=(
      const TrottingConfigurationSpaceCost&) = default;

  // Use defalut move constructor.
  TrottingConfigurationSpaceCost(TrottingConfigurationSpaceCost&&) noexcept = default;

  // Use defalut copy operator.
  TrottingConfigurationSpaceCost& operator=(
      TrottingConfigurationSpaceCost&&) noexcept = default;

  bool useKinematics() const override;

  void set_ref(const double t_start, const double t_period, 
               const Eigen::VectorXd& q_standing, const double step_length,
               const TrottingSwingAngles& swing_angles);

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

private:
  int dimq_, dimv_;
  double t_start_, t_period_, step_length_, v_com_;
  Eigen::VectorXd q_standing_, v_ref_, q_weight_, v_weight_, a_weight_, 
                  qf_weight_, vf_weight_, qi_weight_, vi_weight_, dvi_weight_;
  TrottingSwingAngles swing_angles_;


  void update_q_ref(const double t, Eigen::VectorXd& q_ref) const {
    assert(q_ref.size() == dimq_);
    if (t > t_start_) {
      const double tau = t - t_start_;
      const int steps = std::floor(tau/t_period_);
      const double tau_step = tau - steps * t_period_;
      const double rate = tau_step / t_period_;
      const double sin = std::sin(M_PI*rate);
      const double sin2 = std::sin(M_PI_2*rate);
      if (steps % 2 == 0) {
        q_ref = q_standing_;
        q_ref.coeffRef( 0) += (steps + rate) * step_length_;
        // q_ref.coeffRef( 8) -= swing_angles_.front_swing_thigh;
        q_ref.coeffRef( 9) -= sin2 * swing_angles_.front_swing_knee;
        // q_ref.coeffRef(11) += rate * swing_angles_.hip_stance_thigh;
        // q_ref.coeffRef(12) += rate * swing_angles_.hip_stance_knee;
        // q_ref.coeffRef(14) += rate * swing_angles_.front_stance_thigh;
        // q_ref.coeffRef(15) -= rate * swing_angles_.front_stance_knee;
        // q_ref.coeffRef(17) += rate * swing_angles_.hip_swing_thigh;
        q_ref.coeffRef(18) += sin2 * swing_angles_.hip_swing_knee;
      }
      else {
        q_ref = q_standing_;
        q_ref.coeffRef( 0) += (steps + rate) * step_length_;
        // q_ref.coeffRef( 8) += rate * swing_angles_.front_stance_thigh;
        // q_ref.coeffRef( 9) += rate * swing_angles_.front_stance_knee;
        // q_ref.coeffRef(11) -= rate * swing_angles_.hip_swing_thigh;
        q_ref.coeffRef(12) += sin2 * swing_angles_.hip_swing_knee;
        // q_ref.coeffRef(14) += rate * swing_angles_.front_swing_thigh;
        // q_ref.coeffRef(15) += rate * swing_angles_.front_swing_knee;
        q_ref.coeffRef(15) -= sin2 * swing_angles_.front_swing_knee;
        // q_ref.coeffRef(17) += rate * swing_angles_.hip_stance_thigh;
        // q_ref.coeffRef(18) -= rate * swing_angles_.hip_stance_knee;
      }
    }
    else {
      q_ref = q_standing_;
    }
  }

};

} // namespace idocp

#endif // IDOCP_TROTTING_CONFIGURATION_SPACE_COST_HPP_ 