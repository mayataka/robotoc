#ifndef IDOCP_JOINT_SPACE_COST_HPP_
#define IDOCP_JOINT_SPACE_COST_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {

class JointSpaceCost {
public:
  JointSpaceCost(const Robot& robot, const Eigen::VectorXd& q_weight,  
                 const Eigen::VectorXd& v_wjight, 
                 const Eigen::VectorXd& a_weight,  
                 const Eigen::VectorXd& u_weight,
                 const Eigen::VectorXd& qf_weight,  
                 const Eigen::VectorXd& vf_weight);

  JointSpaceCost(const Robot& robot, const Eigen::VectorXd& q_ref, 
                 const Eigen::VectorXd& q_weight, const Eigen::VectorXd& v_ref, 
                 const Eigen::VectorXd& v_weight, const Eigen::VectorXd& a_ref, 
                 const Eigen::VectorXd& a_weight, const Eigen::VectorXd& u_ref, 
                 const Eigen::VectorXd& u_weight, const Eigen::VectorXd& qf_ref, 
                 const Eigen::VectorXd& qf_weight, 
                 const Eigen::VectorXd& vf_ref, 
                 const Eigen::VectorXd& vf_weight);

  double l(const Robot& robot, const double dtau, const Eigen::VectorXd& q, 
           const Eigen::VectorXd& v, const Eigen::VectorXd& a,
           const Eigen::VectorXd& u);

  double phi(const Robot& robot, const Eigen::VectorXd& q, 
             const Eigen::VectorXd& v);

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

  void phiq(const Robot& robot, const Eigen::VectorXd& q, 
            Eigen::VectorXd& phiq);

  void phiv(const Robot& robot, const Eigen::VectorXd& v, 
            Eigen::VectorXd& phiv);

  void phiqq(const Robot& robot, Eigen::MatrixXd& phiqq);

  void phivv(const Robot& robot, Eigen::MatrixXd& phivv);

  // Prohibits copy constructor.
  JointSpaceCost(const JointSpaceCost&) = delete;

  // Prohibits copy operator.
  JointSpaceCost& operator=(const JointSpaceCost&) = delete;

private:
  unsigned int dimq_, dimv_;
  Eigen::VectorXd q_ref_, v_ref_, a_ref_, u_ref_, q_weight_, v_weight_, 
                  a_weight_, u_weight_;
  Eigen::VectorXd qf_ref_, vf_ref_, qf_weight_, vf_weight_;

};

} // namespace idocp


#endif // IDOCP_JOINT_SPACE_COST_HPP_