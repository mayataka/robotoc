#ifndef IDOCP_TASK_SPACE_6D_COST_HPP_
#define IDOCP_TASK_SPACE_6D_COST_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class TaskSpace6DCost final : public CostFunctionComponentBase {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  TaskSpace6DCost(const Robot& robot, const int end_effector_frame_id);

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

  void set_q_ref(const Eigen::Vector3d& q_ref);

  void set_v_ref(const Eigen::Vector3d& v_ref);

  void set_a_ref(const Eigen::Vector3d& a_ref);
 
  void set_q_weight(const Eigen::Vector3d& q_weight);

  void set_v_weight(const Eigen::Vector3d& v_weight);

  void set_a_weight(const Eigen::Vector3d& a_weight);

  void set_qf_weight(const Eigen::Vector3d& qf_weight);

  void set_vf_weight(const Eigen::Vector3d& vf_weight);

  double l(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s) const override;

  double phi(Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s) const override; 

  void lq(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const override;

  void lv(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const override;

  void la(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s,
          KKTResidual& kkt_residual) const override;

  void lf(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const override {}

  void lqq(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const override;

  void lvv(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const override;

  void laa(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const override;

  void lff(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const override {}

  void phiq(Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, KKTResidual& kkt_residual) const override;

  void phiv(Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, KKTResidual& kkt_residual) const override;

  void phiqq(Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s, KKTMatrix& kkt_matrix) const override;

  void phivv(Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s, KKTMatrix& kkt_matrix) const override;

  void lu(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::VectorXd& u, 
          Eigen::VectorXd& lu) const override {}

  void luu(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::VectorXd& u, 
           Eigen::MatrixXd& Quu) const override {}

private:
  int dimq_, dimv_;
  Eigen::Matrix<double, 6, 1> q_ref_, v_ref_, a_ref_, q_weight_, v_weight_, 
                              a_weight_, qf_weight_, vf_weight_;

};

} // namespace idocp


#endif // IDOCP_TASK_SPACE_6D_COST_HPP_