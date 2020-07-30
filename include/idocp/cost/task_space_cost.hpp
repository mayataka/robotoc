#ifndef IDOCP_TASK_SPACE_COST_HPP_
#define IDOCP_TASK_SPACE_COST_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {

class TaskSpaceCost {
public:
  TaskSpaceCost(const Robot& robot, 
                const Eigen::VectorXd& q_weight,  
                const Eigen::VectorXd& v_weight, 
                const Eigen::VectorXd& a_weight,  
                const Eigen::VectorXd& u_weight,
                const Eigen::VectorXd& qf_weight,  
                const Eigen::VectorXd& vf_weight);

  // Use defalut copy constructor.
  TaskSpaceCost(const TaskSpaceCost&) = default;

  // Use defalut copy operator.
  TaskSpaceCost& operator=(const TaskSpaceCost&) = default;

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

  double l(const Robot& robot, const double dtau, const Eigen::VectorXd& q, 
           const Eigen::VectorXd& v, const Eigen::VectorXd& a,
           const Eigen::VectorXd& u);

  void lq(const Robot& robot, const double dtau, const Eigen::VectorXd& q, 
          Eigen::VectorXd& lq);

  void lv(const Robot& robot, const double dtau, const Eigen::VectorXd& v, 
          Eigen::VectorXd& lv);

  void la(const Robot& robot, const double dtau, const Eigen::VectorXd& a, 
          Eigen::VectorXd& la);

  void lu(const Robot& robot, const double dtau, const Eigen::VectorXd& u, 
          Eigen::VectorXd& lu);

  void lqq(const Robot& robot, const double dtau, Eigen::MatrixXd& lqq);

  void lvv(const Robot& robot, const double dtau, Eigen::MatrixXd& lvv);

  void laa(const Robot& robot, const double dtau, Eigen::MatrixXd& laa);

  void luu(const Robot& robot, const double dtau, Eigen::MatrixXd& luu);

  double phi(const Robot& robot, const Eigen::VectorXd& q, 
             const Eigen::VectorXd& v);

  void phiq(const Robot& robot, const Eigen::VectorXd& q, 
            Eigen::VectorXd& phiq);

  void phiv(const Robot& robot, const Eigen::VectorXd& v, 
            Eigen::VectorXd& phiv);

  void phiqq(const Robot& robot, Eigen::MatrixXd& phiqq);

  void phivv(const Robot& robot, Eigen::MatrixXd& phivv);

private:
  int dimq_, dimv_;
  Eigen::VectorXd q_ref_, v_ref_, a_ref_, u_ref_, q_weight_, v_weight_, 
                  a_weight_, u_weight_, qf_weight_, vf_weight_;

};

} // namespace idocp


#endif // IDOCP_TASK_SPACE_COST_HPP_