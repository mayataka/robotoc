#ifndef IDOCP_JOINT_VELOCITY_UPPER_LIMITS_BARRIER_HPP_
#define IDOCP_JOINT_VELOCITY_UPPER_LIMITS_BARRIER_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {
namespace barrier {

class JointVelocityUpperLimits {
public:
  JointVelocityUpperLimits(const Robot& robot, const double barrier);

  // Use default copy constructor.
  JointVelocityUpperLimits(const JointVelocityUpperLimits&) = default;

  // Use default copy operator.
  JointVelocityUpperLimits& operator=(const JointVelocityUpperLimits&) = default;

  void augmentBarrier(const Robot& robot, const double dtau, 
                      const Eigen::VectorXd& v, Eigen::MatrixXd& Cvv, 
                      Eigen::VectorXd& Cv);

  void augmentBarrierResidual(const Robot& robot, const double dtau, 
                              const Eigen::VectorXd& v, Eigen::VectorXd& Cv);

private:
  unsigned int dimq_, dimv_, dimc_;
  double barrier_;
  Eigen::VectorXd vmax_, residual_;
};

} // namespace barrier
} // namespace idocp


#endif // IDOCP_JOINT_VELOCITY_UPPER_LIMITS_BARRIER_HPP_