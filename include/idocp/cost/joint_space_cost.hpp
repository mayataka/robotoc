#ifndef IDOCP_JOINT_SPACE_COST_HPP_
#define IDOCP_JOINT_SPACE_COST_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

class JointSpaceCost final : public CostFunctionComponentBase {
public:
  using Vector6d = Eigen::Matrix<double, 6, 1>;

  JointSpaceCost(const Robot& robot);

  JointSpaceCost();

  ~JointSpaceCost();

  // Use defalut copy constructor.
  JointSpaceCost(const JointSpaceCost&) = default;

  // Use defalut copy operator.
  JointSpaceCost& operator=(const JointSpaceCost&) = default;

  // Use defalut move constructor.
  JointSpaceCost(JointSpaceCost&&) noexcept = default;

  // Use defalut move assign operator.
  JointSpaceCost& operator=(JointSpaceCost&&) noexcept = default;

  bool useKinematics() const override;

  void set_q_ref(const Eigen::VectorXd& q_ref);

  void set_v_ref(const Eigen::VectorXd& v_ref);

  void set_a_ref(const Eigen::VectorXd& a_ref);
 
  void set_u_ref(const Eigen::VectorXd& u_ref);

  void set_q_weight(const Eigen::VectorXd& q_weight);

  void set_v_weight(const Eigen::VectorXd& v_weight);

  void set_a_weight(const Eigen::VectorXd& a_weight);

  void set_u_weight(const Eigen::VectorXd& u_weight);

  void set_u_passive_weight(const Vector6d& u_passive_weight);

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
          SplitKKTResidual& kkt_residual) const override;

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
           SplitKKTMatrix& kkt_matrix) const override;

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
  int dimq_, dimv_, dimu_;
  Eigen::VectorXd q_ref_, v_ref_, a_ref_, u_ref_, q_weight_, v_weight_, 
                  a_weight_, u_weight_, qf_weight_, vf_weight_;
  Vector6d u_passive_weight_;
};

} // namespace idocp


#endif // IDOCP_JOINT_SPACE_COST_HPP_