#ifndef IDOCP_TASK_SPACE_COST_HPP_
#define IDOCP_TASK_SPACE_COST_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_data.hpp"


namespace idocp {

class TaskSpaceCost {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  TaskSpaceCost(const Robot& robot, const Eigen::VectorXd& q_weight,  
                const Eigen::VectorXd& v_weight, 
                const Eigen::VectorXd& a_weight,  
                const Eigen::VectorXd& u_weight,
                const Eigen::VectorXd& qf_weight,  
                const Eigen::VectorXd& vf_weight);

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

  void set_q_ref(const Eigen::VectorXd& q_ref);

  void set_v_ref(const Eigen::VectorXd& v_ref);

  void set_a_ref(const Eigen::VectorXd& a_ref);
 
  void set_u_ref(const Eigen::VectorXd& u_ref);

  void set_q_weight(const Eigen::VectorXd& q_weight);

  void set_v_weight(const Eigen::VectorXd& v_weight);

  void set_a_weight(const Eigen::VectorXd& a_weight);

  void set_u_weight(const Eigen::VectorXd& u_weight);

  void set_qf_weight(const Eigen::VectorXd& qf_weight);

  void set_vf_weight(const Eigen::VectorXd& vf_weight);

  double l(const Robot& robot, CostFunctionData& data, const double t,
           const double dtau, const Eigen::VectorXd& q, 
           const Eigen::VectorXd& v, const Eigen::VectorXd& a,
           const Eigen::VectorXd& u) const;

  void lq(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::VectorXd& q, 
          Eigen::VectorXd& lq) const;

  void lv(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::VectorXd& v, 
          Eigen::VectorXd& lv) const;

  void la(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::VectorXd& a, 
          Eigen::VectorXd& la) const;

  void lqq(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, Eigen::MatrixXd& lqq) const;

  void lvv(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, Eigen::MatrixXd& lvv) const;

  void laa(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, Eigen::MatrixXd& laa) const;

  double phi(const Robot& robot, CostFunctionData& data, const double t, 
             const Eigen::VectorXd& q, const Eigen::VectorXd& v) const;

  void phiq(const Robot& robot, CostFunctionData& data, const double t, 
            const Eigen::VectorXd& q, Eigen::VectorXd& phiq) const;

  void phiv(const Robot& robot, CostFunctionData& data, const double t, 
            const Eigen::VectorXd& v, Eigen::VectorXd& phiv) const;

  void phiqq(const Robot& robot, CostFunctionData& data, const double t, 
             Eigen::MatrixXd& phiqq) const;

  void phivv(const Robot& robot, CostFunctionData& data, const double t, 
             Eigen::MatrixXd& phivv) const;

private:
  int dimq_, dimv_;
  Eigen::VectorXd q_ref_, v_ref_, a_ref_, u_ref_, q_weight_, v_weight_, 
                  a_weight_, u_weight_, qf_weight_, vf_weight_;

};

} // namespace idocp


#endif // IDOCP_TASK_SPACE_COST_HPP_