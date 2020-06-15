#ifndef IDOCP_JOINT_TORQUE_LOWER_LIMITS_SOFT_HPP_
#define IDOCP_JOINT_TORQUE_LOWER_LIMITS_SOFT_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {
namespace softconstraints {

class JointTorqueLowerLimits {
public:
  JointTorqueLowerLimits(const Robot& robot, const double barrier);

  // Use default copy constructor.
  JointTorqueLowerLimits(const JointTorqueLowerLimits&) = default;

  // Use default copy operator.
  JointTorqueLowerLimits& operator=(const JointTorqueLowerLimits&) = default;

  void setSlackAndDual(const Robot& robot, const double dtau,
                       const Eigen::VectorXd& u);

  void linearizeConstraint(const Robot& robot, const double dtau, 
                           const Eigen::VectorXd& u);

  void updateSlackAndDual(const Robot& robot, const double dtau, 
                          const Eigen::VectorXd& du); 

  void condenseSlackAndDual(const Robot& robot, const double dtau, 
                            Eigen::MatrixXd& Cuu, Eigen::VectorXd& Cu) const;

  void augmentDualResidual(const Robot& robot, const double dtau,
                           Eigen::VectorXd& Cu);

  double squaredConstraintsResidualNrom(const Robot& robot, const double dtau,
                                        const Eigen::VectorXd& u);

private:
  unsigned int dimq_, dimv_, dimc_;
  double barrier_;
  Eigen::VectorXd umin_, slack_, dual_, residual_, FB_residual_, dslack_, 
                  ddual_;
};

} // namespace softconstraints
} // namespace idocp


#endif // IDOCP_JOINT_TORQUE_LOWER_LIMITS_SOFT_HPP_