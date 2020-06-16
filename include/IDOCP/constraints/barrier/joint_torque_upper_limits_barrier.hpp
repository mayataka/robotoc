#ifndef IDOCP_JOINT_TORQUE_UPPER_LIMITS_BARRIER_HPP_
#define IDOCP_JOINT_TORQUE_UPPER_LIMITS_BARRIER_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {
namespace barrier {

class JointTorqueUpperLimits {
public:
  JointTorqueUpperLimits(const Robot& robot, const double barrier);

  // Use default copy constructor.
  JointTorqueUpperLimits(const JointTorqueUpperLimits&) = default;

  // Use default copy operator.
  JointTorqueUpperLimits& operator=(const JointTorqueUpperLimits&) = default;

  void augmentBarrier(const Robot& robot, const double dtau, 
                      const Eigen::VectorXd& u, Eigen::MatrixXd& Cuu, 
                      Eigen::VectorXd& Cu);

  void augmentBarrierResidual(const Robot& robot, const double dtau, 
                              const Eigen::VectorXd& u, Eigen::VectorXd& Cu);

private:
  unsigned int dimq_, dimv_, dimc_;
  double barrier_;
  Eigen::VectorXd umax_, residual_;
};

} // namespace barrier
} // namespace idocp


#endif // IDOCP_JOINT_TORQUE_UPPER_LIMITS_BARRIER_HPP_