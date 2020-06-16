#ifndef IDOCP_JOINT_VELOCITY_UPPER_LIMITS_PDIPM_HPP_
#define IDOCP_JOINT_VELOCITY_UPPER_LIMITS_PDIPM_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {
namespace pdipm {

class JointVelocityUpperLimits {
public:
  JointVelocityUpperLimits(const Robot& robot, const double barrier);

  // Use default copy constructor.
  JointVelocityUpperLimits(const JointVelocityUpperLimits&) = default;

  // Use default copy operator.
  JointVelocityUpperLimits& operator=(const JointVelocityUpperLimits&) = default;

  void setSlackAndDual(const Robot& robot, const double dtau,
                       const Eigen::VectorXd& v);

  void condenseSlackAndDual(const Robot& robot, const double dtau, 
                            const Eigen::VectorXd& v, Eigen::MatrixXd& Cvv, 
                            Eigen::VectorXd& Cv);

  double computeDirectionAndMaxStepSize(const Robot& robot, const double dtau, 
                                        const Eigen::VectorXd& dv);

  void updateSlackAndDual(const Robot& robot, const double step_size);

  void augmentDualResidual(const Robot& robot, const double dtau,
                           Eigen::VectorXd& Cv);

  double squaredConstraintsResidualNrom(const Robot& robot, const double dtau,
                                        const Eigen::VectorXd& v);

private:
  unsigned int dimq_, dimv_, dimc_;
  double barrier_;
  Eigen::VectorXd vmax_, slack_, dual_, residual_, dslack_, ddual_;
};

} // namespace pdipm
} // namespace idocp


#endif // IDOCP_JOINT_VELOCITY_UPPER_LIMITS_PDIPM_HPP_