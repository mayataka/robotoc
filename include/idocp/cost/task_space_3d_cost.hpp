#ifndef IDOCP_TASK_SPACE_3D_COST_HPP_
#define IDOCP_TASK_SPACE_3D_COST_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class TaskSpace3DCost final : public CostFunctionComponentBase {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  TaskSpace3DCost(const Robot& robot, const int frame_id);

  TaskSpace3DCost();

  ~TaskSpace3DCost();

  // Use defalut copy constructor.
  TaskSpace3DCost(const TaskSpace3DCost&) = default;

  // Use defalut copy operator.
  TaskSpace3DCost& operator=(const TaskSpace3DCost&) = default;

  // Use defalut move constructor.
  TaskSpace3DCost(TaskSpace3DCost&&) noexcept = default;

  // Use defalut copy operator.
  TaskSpace3DCost& operator=(TaskSpace3DCost&&) noexcept = default;

  bool useKinematics() const override;

  void set_q_3d_ref(const Eigen::Vector3d& q_3d_ref);

  void set_q_3d_weight(const Eigen::Vector3d& q_3d_weight);

  void set_qf_3d_weight(const Eigen::Vector3d& qf_3d_weight);

  double l(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s) const override;

  double phi(Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s) const override; 

  void lq(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const override;

  void lv(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const override {}

  void la(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s,
          KKTResidual& kkt_residual) const override {}

  void lf(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const override {}

  void lqq(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const override;

  void lvv(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const override {}

  void laa(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const override {}

  void lff(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const override {}

  void phiq(Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, KKTResidual& kkt_residual) const override;

  void phiv(Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, KKTResidual& kkt_residual) const override {}

  void phiqq(Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s, KKTMatrix& kkt_matrix) const override;

  void phivv(Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s, KKTMatrix& kkt_matrix) const override {}

  void lu(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::VectorXd& u, 
          Eigen::VectorXd& lu) const override {}

  void luu(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::VectorXd& u, 
           Eigen::MatrixXd& Quu) const override {}

private:
  int frame_id_;
  Eigen::Vector3d q_3d_ref_, q_3d_weight_, qf_3d_weight_;

};

} // namespace idocp


#endif // IDOCP_TASK_SPACE_3D_COST_HPP_