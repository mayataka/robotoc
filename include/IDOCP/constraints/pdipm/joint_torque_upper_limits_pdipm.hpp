#ifndef IDOCP_JOINT_TORQUE_UPPER_LIMITS_PDIPM_HPP_
#define IDOCP_JOINT_TORQUE_UPPER_LIMITS_PDIPM_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {
namespace pdipm {

class JointTorqueUpperLimits {
public:
  JointTorqueUpperLimits(const Robot& robot, const double barrier);

  // Use default copy constructor.
  JointTorqueUpperLimits(const JointTorqueUpperLimits&) = default;

  // Use default copy operator.
  JointTorqueUpperLimits& operator=(const JointTorqueUpperLimits&) = default;

  void setSlackAndDual(const Robot& robot, const double dtau,
                       const Eigen::VectorXd& u);

  void condenseSlackAndDual(const Robot& robot, const double dtau, 
                            const Eigen::VectorXd& u, Eigen::MatrixXd& Cuu, 
                            Eigen::VectorXd& Cu);

  double computeDirectionAndMaxStepSize(const Robot& robot, const double dtau, 
                                        const Eigen::VectorXd& du);

  void updateSlackAndDual(const Robot& robot, const double step_size);

  void augmentDualResidual(const Robot& robot, const double dtau,
                           Eigen::VectorXd& Cu);

  double squaredConstraintsResidualNrom(const Robot& robot, const double dtau,
                                        const Eigen::VectorXd& u);

private:
  unsigned int dimq_, dimv_, dimc_;
  double barrier_;
  Eigen::VectorXd umax_, slack_, dual_, residual_, dslack_, ddual_;
};

} // namespace pdipm
} // namespace idocp


#endif // IDOCP_JOINT_TORQUE_UPPER_LIMITS_PDIPM_HPP_