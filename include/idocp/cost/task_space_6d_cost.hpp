#ifndef IDOCP_TASK_SPACE_6D_COST_HPP_
#define IDOCP_TASK_SPACE_6D_COST_HPP_

#include "Eigen/Core"
#include "pinocchio/spatial/se3.hpp"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

class TaskSpace6DCost final : public CostFunctionComponentBase {
public:
  using Vector6d = Eigen::Matrix<double, 6, 1>;

  TaskSpace6DCost(const Robot& robot, const int frame_id);

  TaskSpace6DCost();

  ~TaskSpace6DCost();

  // Use defalut copy constructor.
  TaskSpace6DCost(const TaskSpace6DCost&) = default;

  // Use defalut copy operator.
  TaskSpace6DCost& operator=(const TaskSpace6DCost&) = default;

  // Use defalut move constructor.
  TaskSpace6DCost(TaskSpace6DCost&&) noexcept = default;

  // Use defalut copy operator.
  TaskSpace6DCost& operator=(TaskSpace6DCost&&) noexcept = default;

  bool useKinematics() const override;

  void set_q_6d_ref(const pinocchio::SE3& SE3_ref);

  void set_q_6d_ref(const Eigen::Vector3d& position_ref, 
                    const Eigen::Matrix3d& rotation_mat_ref);

  void set_q_6d_weight(const Eigen::Vector3d& position_weight, 
                       const Eigen::Vector3d& rotation_weight);

  void set_q_6d_weight(const Vector6d& q_6d_weight);

  void set_qf_6d_weight(const Eigen::Vector3d& position_weight, 
                        const Eigen::Vector3d& rotation_weight);

  void set_qf_6d_weight(const Vector6d& q_6d_weight);

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
  pinocchio::SE3 SE3_ref_, SE3_ref_inv_;
  Eigen::VectorXd q_6d_weight_, qf_6d_weight_;

};

} // namespace idocp


#endif // IDOCP_TASK_SPACE_6D_COST_HPP_