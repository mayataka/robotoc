#ifndef IDOCP_JOINT_SPACE_SOFT_CONSTRAINTS_HPP_
#define IDOCP_JOINT_SPACE_SOFT_CONSTRAINTS_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"
#include "constraints/soft/joint_position_upper_limits_soft.hpp"
#include "constraints/soft/joint_position_lower_limits_soft.hpp"
#include "constraints/soft/joint_velocity_upper_limits_soft.hpp"
#include "constraints/soft/joint_velocity_lower_limits_soft.hpp"
#include "constraints/soft/joint_torque_upper_limits_soft.hpp"
#include "constraints/soft/joint_torque_lower_limits_soft.hpp"


namespace idocp {

class JointSpaceSoftConstraints {
public:
  JointSpaceSoftConstraints(const Robot& robot);

  // Use default copy constructor.
  JointSpaceSoftConstraints(const JointSpaceSoftConstraints&) = default;

  // Use default copy operator.
  JointSpaceSoftConstraints& operator=(const JointSpaceSoftConstraints&) = default;

  void setSlackAndDual(const Robot& robot, const double dtau, 
                       const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                       const Eigen::VectorXd& u, const Eigen::VectorXd& a);

  void linearizeConstraint(const Robot& robot, const double dtau,
                           const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                           const Eigen::VectorXd& a, const Eigen::VectorXd& u);

  void updateSlackAndDual(const Robot& robot, const double dtau, 
                          const Eigen::VectorXd& dq, const Eigen::VectorXd& dv, 
                          const Eigen::VectorXd& da, const Eigen::VectorXd& du);

  void condenseSlackAndDual(const Robot& robot, const double dtau, 
                            Eigen::MatrixXd& Cqq,  Eigen::MatrixXd& Cvv, 
                            Eigen::MatrixXd& Caa, Eigen::MatrixXd& Cuu, 
                            Eigen::VectorXd& Cq, Eigen::VectorXd& Cv, 
                            Eigen::VectorXd& Ca, Eigen::VectorXd& Cu);

  void condenseSlackAndDual(const Robot& robot, const double dtau, 
                            Eigen::MatrixXd& Cqq,  Eigen::MatrixXd& Cvv, 
                            Eigen::MatrixXd& Caa,  Eigen::VectorXd& Cq, 
                            Eigen::VectorXd& Cv, Eigen::VectorXd& Ca);

  void condenseSlackAndDual(const Robot& robot, const double dtau, 
                            Eigen::MatrixXd& Cuu, Eigen::VectorXd& Cu);

  void augmentDualResidual(const Robot& robot, const double dtau, 
                           Eigen::VectorXd& Cq,
                           Eigen::VectorXd& Cv, 
                           Eigen::VectorXd& Ca, 
                           Eigen::VectorXd& Cu);

  double squaredConstraintsResidualNrom(const Robot& robot, const double dtau,
                                        const Eigen::VectorXd& q, 
                                        const Eigen::VectorXd& v, 
                                        const Eigen::VectorXd& a, 
                                        const Eigen::VectorXd& u);

private:
  double barrier_;
  softconstraints::JointPositionUpperLimits position_upper_limits_;
  softconstraints::JointPositionLowerLimits position_lower_limits_;
  softconstraints::JointVelocityUpperLimits velocity_upper_limits_;
  softconstraints::JointVelocityLowerLimits velocity_lower_limits_;
  softconstraints::JointTorqueUpperLimits torque_upper_limits_;
  softconstraints::JointTorqueLowerLimits torque_lower_limits_;
};

} // namespace idocp


#endif // IDOCP_JOINT_SPACE_CONSTRAINTS_HPP_