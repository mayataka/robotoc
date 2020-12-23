#ifndef IDOCP_TROTTING_CONFIGURATION_COST_HPP_
#define IDOCP_TROTTING_CONFIGURATION_COST_HPP_

#include <cmath>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

class TrottingConfigurationCost final : public CostFunctionComponentBase {
public:

  TrottingConfigurationCost(const Robot& robot);

  TrottingConfigurationCost();

  ~TrottingConfigurationCost();

  // Use defalut copy constructor.
  TrottingConfigurationCost(const TrottingConfigurationCost&) = default;

  // Use defalut copy operator.
  TrottingConfigurationCost& operator=(
      const TrottingConfigurationCost&) = default;

  // Use defalut move constructor.
  TrottingConfigurationCost(TrottingConfigurationCost&&) noexcept = default;

  // Use defalut copy operator.
  TrottingConfigurationCost& operator=(
      TrottingConfigurationCost&&) noexcept = default;

  bool useKinematics() const override;

  void set_ref(const double t_start, const double t_period, 
               const Eigen::VectorXd q_standing, const double step_length,
               const double front_thigh_swing_angle=0.3, 
               const double front_knee_swing_angle=0.5,
               const double hip_thigh_swing_angle=0.3, 
               const double hip_knee_swing_angle=0.5);

  void set_q_weight(const Eigen::VectorXd& q_weight);

  void set_v_weight(const Eigen::VectorXd& v_weight);

  void set_a_weight(const Eigen::VectorXd& a_weight);

  void set_qf_weight(const Eigen::VectorXd& qf_weight);

  void set_vf_weight(const Eigen::VectorXd& vf_weight);

  double l(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s) const override;

  double phi(Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s) const override; 

  void lq(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          SplitKKTResidual& kkt_residual) const override;

  void lv(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          SplitKKTResidual& kkt_residual) const override; 

  void la(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s,
          SplitKKTResidual& kkt_residual) const override;

  void lf(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          SplitKKTResidual& kkt_residual) const override {}

  void lu(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          SplitKKTResidual& kkt_residual) const override {}

  void lqq(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           SplitKKTMatrix& kkt_matrix) const override;

  void lvv(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           SplitKKTMatrix& kkt_matrix) const override;

  void laa(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           SplitKKTMatrix& kkt_matrix) const override;

  void lff(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           SplitKKTMatrix& kkt_matrix) const override {}

  void luu(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           SplitKKTMatrix& kkt_matrix) const override {}

  void phiq(Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, 
            SplitKKTResidual& kkt_residual) const override;

  void phiv(Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, 
            SplitKKTResidual& kkt_residual) const override;

  void phiqq(Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s, 
             SplitKKTMatrix& kkt_matrix) const override;

  void phivv(Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s, 
             SplitKKTMatrix& kkt_matrix) const override;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int dimq_, dimv_;
  double t_start_, t_period_, step_length_, v_com_;
  Eigen::VectorXd q_standing_, q_even_step_, q_odd_step_, v_ref_,
                  q_weight_, v_weight_, a_weight_, qf_weight_, vf_weight_;


  void update_q_ref(const double t, Eigen::VectorXd& q_ref) const {
    assert(q_ref.size() == dimq_);
    if (t > t_start_) {
      const int steps = std::floor((t-t_start_)/t_period_);
      if (steps % 2 == 0) {
        q_ref = q_even_step_;
        q_ref.coeffRef(0) += steps * step_length_;
      }
      else {
        q_ref = q_odd_step_;
        q_ref.coeffRef(0) += steps * step_length_;
      }
    }
    else {
      q_ref = q_standing_;
    }
  }

};

} // namespace idocp

#endif // IDOCP_TROTTING_CONFIGURATION_COST_HPP_ 