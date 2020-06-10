#ifndef IDOCP_JOINT_SPACE_CONSTRAINTS_HPP_
#define IDOCP_JOINT_SPACE_CONSTRAINTS_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {

class JointSpaceConstraints {
public:
  JointSpaceConstraints(const Robot& robot);

  void Cq(const Robot& robot, const double dtau, const Eigen::VectorXd& q, 
          Eigen::VectorXd& Cq);

  void Cv(const Robot& robot, const double dtau, const Eigen::VectorXd& v, 
          Eigen::VectorXd& Cv);

  void Cu(const Robot& robot, const double dtau, const Eigen::VectorXd& u, 
          Eigen::VectorXd& Cu);

  // Prohibits copy constructor.
  JointSpaceConstraints(const JointSpaceConstraints&) = delete;

  // Prohibits copy operator.
  JointSpaceConstraints& operator=(const JointSpaceConstraints&) = delete;

private:
  Eigen::VectorXd q_max_, q_min_, v_max_, v_min_, u_max_, u_min_;
  unsigned int dimq_, dimv_;

};

} // namespace idocp

#endif // IDOCP_JOINT_SPACE_CONSTRAINTS_HPP_