#ifndef IDOCP_TROTTING_FOOT_STEP_COST_HPP_
#define IDOCP_TROTTING_FOOT_STEP_COST_HPP_

#include <cmath>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

class TrottingFootStepCost final : public CostFunctionComponentBase {
public:
  TrottingFootStepCost(const Robot& robot, const int foot_frame_id);

  TrottingFootStepCost();

  ~TrottingFootStepCost();

  // Use defalut copy constructor.
  TrottingFootStepCost(const TrottingFootStepCost&) = default;

  // Use defalut copy operator.
  TrottingFootStepCost& operator=(const TrottingFootStepCost&) = default;

  // Use defalut move constructor.
  TrottingFootStepCost(TrottingFootStepCost&&) noexcept = default;

  // Use defalut copy operator.
  TrottingFootStepCost& operator=(TrottingFootStepCost&&) noexcept = default;

  bool useKinematics() const override;

  void set_ref(Robot& robot, const Eigen::VectorXd& q_standing, 
               const double step_length, const double step_height);

  void set_period(const double t_start, const double t_period, 
                  const bool is_initial_foot_step);

  void set_q_3d_weight(const Eigen::Vector3d& q_3d_weight);

  void set_qf_3d_weight(const Eigen::Vector3d& qf_3d_weight);

  void set_qi_3d_weight(const Eigen::Vector3d& qi_3d_weight);

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


  Eigen::Vector3d q_3d_ref(const double t) const {
    if (is_initial_foot_step_) {
      if (t > t_start_+t_period_) {
        const double tau = t - t_start_;
        const int steps = std::floor(tau/t_period_);
        Eigen::Vector3d q_3d_ref_vec(q_3d_ref_init_);
        if (steps % 2 == 0) {
          const double tau_step = tau - steps * t_period_;
          const double rate = tau_step / t_period_;
          q_3d_ref_vec.coeffRef(0) += (0.5*(steps-1) + rate) * step_length_;
          q_3d_ref_vec.coeffRef(2) += std::sin(M_PI*rate) * step_height_;
        }
        else {
          q_3d_ref_vec.coeffRef(0) += 0.5*steps * step_length_;
        }
        return q_3d_ref_vec;
      }
      else if (t > t_start_) {
        Eigen::Vector3d q_3d_ref_vec(q_3d_ref_init_);
        const double tau_step = t - t_start_;
        const double rate = tau_step / t_period_;
        q_3d_ref_vec.coeffRef(0) += 0.5 * rate * step_length_;
        q_3d_ref_vec.coeffRef(2) += std::sin(M_PI*rate) * step_height_;
        return q_3d_ref_vec;
      }
      else {
        return q_3d_ref_init_;
      }
    }
    else {
      if (t > t_start_+t_period_) {
        const double tau = t - t_start_ - t_period_;
        const int steps = std::floor(tau/t_period_);
        Eigen::Vector3d q_3d_ref_vec(q_3d_ref_init_);
        if (steps % 2 == 0) {
          const double tau_step = tau - steps * t_period_;
          const double rate = tau_step / t_period_;
          q_3d_ref_vec.coeffRef(0) += (0.5*steps + rate) * step_length_;
          q_3d_ref_vec.coeffRef(2) += std::sin(M_PI*rate) * step_height_;
        }
        else {
          q_3d_ref_vec.coeffRef(0) += 0.5*(steps+1) * step_length_;
        }
        return q_3d_ref_vec;
      }
      else {
        return q_3d_ref_init_;
      }
    }
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int frame_id_;
  double t_start_, t_period_, step_length_, step_height_;
  bool is_initial_foot_step_;
  Eigen::Vector3d q_3d_ref_init_, q_3d_weight_, qf_3d_weight_, qi_3d_weight_;

};

} // namespace idocp

#endif // IDOCP_TROTTING_FOOT_STEP_COST_HPP_ 