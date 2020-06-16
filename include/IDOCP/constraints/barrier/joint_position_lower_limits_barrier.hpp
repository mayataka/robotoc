#ifndef IDOCP_JOINT_POSITION_LOWER_LIMITS_BARRIER_HPP_
#define IDOCP_JOINT_POSITION_LOWER_LIMITS_BARRIER_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {
namespace barrier {

class JointPositionLowerLimits {
public:
  JointPositionLowerLimits(const Robot& robot, const double barrier);

  // Use default copy constructor.
  JointPositionLowerLimits(const JointPositionLowerLimits&) = default;

  // Use default copy operator.
  JointPositionLowerLimits& operator=(const JointPositionLowerLimits&) = default;

  void augmentBarrier(const Robot& robot, const double dtau, 
                      const Eigen::VectorXd& q, Eigen::MatrixXd& Cqq, 
                      Eigen::VectorXd& Cq);

  void augmentBarrierResidual(const Robot& robot, const double dtau, 
                              const Eigen::VectorXd& q, Eigen::VectorXd& Cq);

private:
  unsigned int dimq_, dimv_, dimc_;
  double barrier_;
  Eigen::VectorXd qmin_, residual_;
};

} // namespace barrier
} // namespace idocp


#endif // IDOCP_JOINT_POSITION_LOWER_LIMITS_BARRIER_HPP_