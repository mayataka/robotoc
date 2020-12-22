#ifndef IDOCP_FOOT_STEP_TROTTING_COST_HPP_
#define IDOCP_FOOT_STEP_TROTTING_COST_HPP_

#include <cmath>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

class FootStepTrottingCost final : public CostFunctionComponentBase {
public:

  FootStepTrottingCost(const Robot& robot, const int foot_frame_id);

  FootStepTrottingCost();

  ~FootStepTrottingCost();

  // Use defalut copy constructor.
  FootStepTrottingCost(const FootStepTrottingCost&) = default;

  // Use defalut copy operator.
  FootStepTrottingCost& operator=(const FootStepTrottingCost&) = default;

  // Use defalut move constructor.
  FootStepTrottingCost(FootStepTrottingCost&&) noexcept = default;

  // Use defalut copy operator.
  FootStepTrottingCost& operator=(FootStepTrottingCost&&) noexcept = default;

  bool useKinematics() const override;

  void set_q_3d_ref(const Eigen::Vector3d& q_3d_ref_init, 
                    const double step_length, const double step_height);

  void set_period(const double t_start, const double t_period);

  void set_q_3d_weight(const Eigen::Vector3d& q_3d_weight);

  void set_qf_3d_weight(const Eigen::Vector3d& qf_3d_weight);

  double l(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s) const override;

  double phi(Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s) const override; 

  void lq(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          SplitKKTResidual& kkt_residual) const override;

  void lv(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          SplitKKTResidual& kkt_residual) const override {}

  void la(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s,
          SplitKKTResidual& kkt_residual) const override {}

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
           SplitKKTMatrix& kkt_matrix) const override {}

  void laa(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           SplitKKTMatrix& kkt_matrix) const override {}

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
            SplitKKTResidual& kkt_residual) const override {}

  void phiqq(Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s, 
             SplitKKTMatrix& kkt_matrix) const override;

  void phivv(Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s, 
             SplitKKTMatrix& kkt_matrix) const override {}

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int frame_id_;
  double t_start_, t_period_, step_length_, step_height_;
  Eigen::Vector3d q_3d_ref_init_, q_3d_weight_, qf_3d_weight_;


  Eigen::Vector3d q_3d_ref(const double t) const {
    if (t > t_start_) {
      const int phase = std::floor((t-t_start_)/t_period_) + 1;
      Eigen::Vector3d q_3d_ref_vec(q_3d_ref_init_);
      if (phase % 2 == 0) {
        q_3d_ref_vec.coeffRef(0) += (phase/2) * step_length_;
      }
      else {
        q_3d_ref_vec.coeffRef(0) += ((phase-1)/2 + 1) * step_length_;
        q_3d_ref_vec.coeffRef(2) += step_height_;
      }
      return q_3d_ref_vec;
    }
    else {
      return q_3d_ref_init_;
    }
  }

};

} // namespace idocp

#endif // IDOCP_FOOT_STEP_TROTTING_COST_HPP_ 