#ifndef IDOCP_JOINT_VELOCITY_LOWER_LIMITS_BARRIER_HPP_
#define IDOCP_JOINT_VELOCITY_LOWER_LIMITS_BARRIER_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {
namespace barrier {

class JointVelocityLowerLimits {
public:
  JointVelocityLowerLimits(const Robot& robot, const double barrier);

  // Use default copy constructor.
  JointVelocityLowerLimits(const JointVelocityLowerLimits&) = default;

  // Use default copy operator.
  JointVelocityLowerLimits& operator=(const JointVelocityLowerLimits&) = default;

  void augmentBarrier(const Robot& robot, const double dtau, 
                      const Eigen::VectorXd& v, Eigen::MatrixXd& Cvv, 
                      Eigen::VectorXd& Cv);

  void augmentBarrierResidual(const Robot& robot, const double dtau, 
                              const Eigen::VectorXd& v, Eigen::VectorXd& Cv);

private:
  unsigned int dimq_, dimv_, dimc_;
  double barrier_;
  Eigen::VectorXd vmin_, residual_;
};

} // namespace barrier
} // namespace idocp


#endif // IDOCP_JOINT_VELOCITY_LOWER_LIMITS_BARRIER_HPP_