#ifndef IDOCP_TROTTING_FOOT_HEIGHT_COST_HPP_
#define IDOCP_TROTTING_FOOT_HEIGHT_COST_HPP_

#include <cmath>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"

namespace idocp {

class TrottingFootHeightCost final : public CostFunctionComponentBase {
public:
  TrottingFootHeightCost(const Robot &robot, const int knee_joint_id, 
                         const int hip_joint_id, const double foot_length, 
                         const double thigh_length);

  TrottingFootHeightCost();

  ~TrottingFootHeightCost();

  // Use defalut copy constructor.
  TrottingFootHeightCost(const TrottingFootHeightCost&) = default;

  // Use defalut copy operator.
  TrottingFootHeightCost&operator=(const TrottingFootHeightCost&) = default;

  // Use defalut move constructor.
  TrottingFootHeightCost(TrottingFootHeightCost&&) noexcept = default;

  // Use defalut copy operator.
  TrottingFootHeightCost&operator=(TrottingFootHeightCost&&) noexcept = default;

  bool useKinematics() const override;

  void set_ref(const Eigen::VectorXd& q_standing, const double swing_height);

  void set_period(const double t_start, const double t_period);

  void set_q_weight(const double q_weight);

  void set_qf_weight(const double qf_weight);

  void set_qi_weight(const double qi_weight);

  double computeStageCost(Robot &robot, CostFunctionData &data, const double t,
                          const double dtau, const SplitSolution &s) const;

  double computeTerminalCost(Robot &robot, CostFunctionData &data,
                              const double t, const SplitSolution &s) const;

  double computeImpulseCost(Robot &robot, CostFunctionData &data,
                            const double t,
                            const ImpulseSplitSolution &s) const;

  void computeStageCostDerivatives(Robot &robot, CostFunctionData &data,
                                    const double t, const double dtau,
                                    const SplitSolution &s,
                                    SplitKKTResidual &kkt_residual) const;

  void computeTerminalCostDerivatives(Robot &robot, CostFunctionData &data,
                                      const double t, const SplitSolution &s,
                                      SplitKKTResidual &kkt_residual) const;

  void computeImpulseCostDerivatives(Robot &robot, CostFunctionData &data,
                                      const double t,
                                      const ImpulseSplitSolution &s,
                                      ImpulseSplitKKTResidual &kkt_residual) const;

  void computeStageCostHessian(Robot &robot, CostFunctionData &data,
                                const double t, const double dtau,
                                const SplitSolution &s,
                                SplitKKTMatrix &kkt_matrix) const;

  void computeTerminalCostHessian(Robot &robot, CostFunctionData &data,
                                  const double t, const SplitSolution &s,
                                  SplitKKTMatrix &kkt_matrix) const;

  void computeImpulseCostHessian(Robot &robot, CostFunctionData &data,
                                  const double t, const ImpulseSplitSolution &s,
                                  ImpulseSplitKKTMatrix &kkt_matrix) const;


  double hdiff(const double t, const Eigen::VectorXd& q) const {
    if (t > t_start_) {
      const int steps = std::floor((t-t_start_)/t_period_);
      if (steps % 2 == 0) 
        return (dist(q) + height_ref(t) - stance_height_);
      else 
        return (dist(q) - stance_height_);
    }
    else {
      return (dist(q) - stance_height_);
    }
  }

  double ddist_dth1(const double t, const Eigen::VectorXd& q) const {
    return (- thigh_length_ * std::sin(th1(q)) - foot_length_ * std::sin(th1(q)+th2(q)));
  }

  double ddist_dth2(const double t, const Eigen::VectorXd& q) const {
    return (- foot_length_ * std::sin(th1(q)+th2(q)));
  }

  double dddist_dth11(const double t, const Eigen::VectorXd& q) const {
    return (- thigh_length_ * std::cos(th1(q)) - foot_length_ * std::cos(th1(q)+th2(q)));
  }

  double dddist_dth12(const double t, const Eigen::VectorXd& q) const {
    return dddist_dth22(t, q);
  }

  double dddist_dth22(const double t, const Eigen::VectorXd& q) const {
    return (- foot_length_ * std::cos(th1(q)+th2(q)));
  }

  double dist(const Eigen::VectorXd& q) const {
    return (thigh_length_ * std::cos(th1(q)) + foot_length_ * std::cos(th1(q)+th2(q)));
  }

  double th1(const Eigen::VectorXd& q) const {
    return q.coeff(hip_joint_id_+1);
  }

  double th2(const Eigen::VectorXd& q) const {
    return q.coeff(knee_joint_id_+1);
  }

  double height_ref(const double t) const {
    const int steps = std::floor((t-t_start_)/t_period_);
    const double tau = t - t_start_ - steps * t_period_;
    return swing_height_ * std::sin(M_PI*tau/t_period_);
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int foot_frame_id_, thigh_frame_id_, knee_joint_id_, hip_joint_id_;
  double foot_length_, thigh_length_, stance_height_, swing_height_, 
         t_start_, t_period_, q_weight_, qf_weight_, qi_weight_;

};

} // namespace idocp

#endif // IDOCP_TROTTING_FOOT_HEIGHT_COST_HPP_