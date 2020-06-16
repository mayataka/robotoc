#ifndef IDOCP_JOINT_POSITION_UPPER_LIMITS_BARRIER_HPP_
#define IDOCP_JOINT_POSITION_UPPER_LIMITS_BARRIER_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {
namespace barrier {

class JointPositionUpperLimits {
public:
  JointPositionUpperLimits(const Robot& robot, const double barrier);

  // Use default copy constructor.
  JointPositionUpperLimits(const JointPositionUpperLimits&) = default;

  // Use default copy operator.
  JointPositionUpperLimits& operator=(const JointPositionUpperLimits&) = default;

  void augmentBarrier(const Robot& robot, const double dtau, 
                      const Eigen::VectorXd& q, Eigen::MatrixXd& Cqq, 
                      Eigen::VectorXd& Cq);

  void augmentBarrierResidual(const Robot& robot, const double dtau, 
                              const Eigen::VectorXd& q, Eigen::VectorXd& Cq);

private:
  unsigned int dimq_, dimv_, dimc_;
  double barrier_;
  Eigen::VectorXd qmax_, residual_;
};

} // namespace barrier
} // namespace idocp


#endif // IDOCP_JOINT_POSITION_UPPER_LIMITS_BARRIER_HPP_