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

  TaskSpaceCost(const Robot& robot, const int end_effector_frame_id);

  TaskSpaceCost();

  ~TaskSpaceCost();

  // Use defalut copy constructor.
  TaskSpaceCost(const TaskSpaceCost&) = default;

  // Use defalut copy operator.
  TaskSpaceCost& operator=(const TaskSpaceCost&) = default;

  // Use defalut move constructor.
  TaskSpaceCost(TaskSpaceCost&&) noexcept = default;

  // Use defalut copy operator.
  TaskSpaceCost& operator=(TaskSpaceCost&&) noexcept = default;

  void set_q_ref(const Eigen::Vector3d& q_ref);

  void set_v_ref(const Eigen::Vector3d& v_ref);

  void set_a_ref(const Eigen::Vector3d& a_ref);
 
  void set_q_weight(const Eigen::Vector3d& q_weight);

  void set_v_weight(const Eigen::Vector3d& v_weight);

  void set_a_weight(const Eigen::Vector3d& a_weight);

  void set_qf_weight(const Eigen::Vector3d& qf_weight);

  void set_vf_weight(const Eigen::Vector3d& vf_weight);

  double l(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s) const override;

  double phi(const Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s) const override; 

  void lq(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const override;

  void lv(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const override;

  void la(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s,
          KKTResidual& kkt_residual) const override;

  void lf(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const override {}

  void lqq(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const override;

  void lvv(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const override;

  void laa(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const override;

  void lff(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const override {}

  void phiq(const Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, KKTResidual& kkt_residual) const override;

  void phiv(const Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, KKTResidual& kkt_residual) const override;

  void phiqq(const Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s, KKTMatrix& kkt_matrix) const override;

  void phivv(const Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s, KKTMatrix& kkt_matrix) const override;

  void lu(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::VectorXd& u, 
          Eigen::VectorXd& lu) const override {}

  void luu(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::VectorXd& u, 
           Eigen::MatrixXd& Quu) const override {}

private:
  int dimq_, dimv_;
  Eigen::Vector3d q_ref_, v_ref_, a_ref_, q_weight_, v_weight_, a_weight_, 
                  qf_weight_, vf_weight_;

};

} // namespace idocp


#endif // IDOCP_TASK_SPACE_3D_COST_HPP_