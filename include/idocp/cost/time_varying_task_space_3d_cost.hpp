#ifndef IDOCP_TIME_VARYING_TASK_SPACE_3D_COST_HPP_
#define IDOCP_TIME_VARYING_TASK_SPACE_3D_COST_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

class TimeVaryingTaskSpace3DCost final : public CostFunctionComponentBase {
public:
  TimeVaryingTaskSpace3DCost(const Robot& robot, const int frame_id);

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

  void set_ref(const double t0, const Eigen::Vector3d& q0, 
               const Eigen::Vector3d& v0);

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
  double t0_;
  Eigen::Vector3d q0_, v0_, q_3d_weight_, qf_3d_weight_;

};

} // namespace idocp


#endif // IDOCP_TIME_VARYING_TASK_SPACE_3D_COST_HPP_ 