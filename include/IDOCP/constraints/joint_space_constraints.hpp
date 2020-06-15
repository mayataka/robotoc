#ifndef IDOCP_JOINT_SPACE_CONSTRAINTS_HPP_
#define IDOCP_JOINT_SPACE_CONSTRAINTS_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"
#include "constraints/joint_position_upper_limits.hpp"
#include "constraints/joint_position_lower_limits.hpp"
#include "constraints/joint_velocity_upper_limits.hpp"
#include "constraints/joint_velocity_lower_limits.hpp"
#include "constraints/joint_torque_upper_limits.hpp"
#include "constraints/joint_torque_lower_limits.hpp"


namespace idocp {

class JointSpaceConstraints {
public:
  JointSpaceConstraints(const Robot& robot);

  // Use default copy constructor.
  JointSpaceConstraints(const JointSpaceConstraints&) = default;

  // Use default copy operator.
  JointSpaceConstraints& operator=(const JointSpaceConstraints&) = default;

  void setSlackAndDual(const Robot& robot, const double dtau, 
                       const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                       const Eigen::VectorXd& u);

  void linearizeConstraint(const Robot& robot, const double dtau,
                           const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                           const Eigen::VectorXd& u);

  void updateSlackAndDual(const Robot& robot, const double dtau, 
                          const Eigen::MatrixXd& du_dq,
                          const Eigen::MatrixXd& du_dv,
                          const Eigen::MatrixXd& du_da,
                          const Eigen::VectorXd& dq, const Eigen::VectorXd& dv, 
                          const Eigen::VectorXd& da);

  void condenseSlackAndDual(const Robot& robot, const double dtau, 
                            const Eigen::MatrixXd& du_dq, 
                            const Eigen::MatrixXd& du_dv,
                            const Eigen::MatrixXd& du_da, Eigen::MatrixXd& Cqq, 
                            Eigen::MatrixXd& Cqv, Eigen::MatrixXd& Cqa, 
                            Eigen::MatrixXd& Cvv, Eigen::MatrixXd& Cva, 
                            Eigen::MatrixXd& Caa, Eigen::VectorXd& Cq, 
                            Eigen::VectorXd& Cv, Eigen::VectorXd& Ca);

  void augmentDualResidual(const Robot& robot, const double dtau, 
                           const Eigen::MatrixXd& du_dq,
                           const Eigen::MatrixXd& du_dv, 
                           const Eigen::MatrixXd& du_da, Eigen::VectorXd& Cq,
                           Eigen::VectorXd& Cv, Eigen::VectorXd& Ca);

  double squaredConstraintsResidualNrom(const Robot& robot, const double dtau,
                                        const Eigen::VectorXd& q, 
                                        const Eigen::VectorXd& v, 
                                        const Eigen::VectorXd& u);

private:
  double barrier_;
  JointPositionUpperLimits position_upper_limits_;
  JointPositionLowerLimits position_lower_limits_;
  JointVelocityUpperLimits velocity_upper_limits_;
  JointVelocityLowerLimits velocity_lower_limits_;
  JointTorqueUpperLimits torque_upper_limits_;
  JointTorqueLowerLimits torque_lower_limits_;
};

} // namespace idocp


#endif // IDOCP_JOINT_SPACE_CONSTRAINTS_HPP_