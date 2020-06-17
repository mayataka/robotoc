#ifndef IDOCP_JOINT_TORQUE_LOWER_LIMITS_PDIPM_HPP_
#define IDOCP_JOINT_TORQUE_LOWER_LIMITS_PDIPM_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {
namespace pdipm {

class JointTorqueLowerLimits {
public:
  JointTorqueLowerLimits(const Robot& robot, const double barrier);

  // Use default copy constructor.
  JointTorqueLowerLimits(const JointTorqueLowerLimits&) = default;

  // Use default copy operator.
  JointTorqueLowerLimits& operator=(const JointTorqueLowerLimits&) = default;

  bool isFeasible(const Robot& robot, const Eigen::VectorXd& u);

  void setSlackAndDual(const Robot& robot, const double dtau,
                       const Eigen::VectorXd& u);

  void condenseSlackAndDual(const Robot& robot, const double dtau, 
                            const Eigen::VectorXd& u, Eigen::MatrixXd& Cuu, 
                            Eigen::VectorXd& Cu);

  std::pair<double, double> computeDirectionAndMaxStepSize(
      const Robot& robot, const double fraction_to_boundary_rate, 
      const double dtau, const Eigen::VectorXd& dq);

  void updateSlack(const double step_size);

  void updateDual(const double step_size);

  double slackBarrier();

  double slackBarrier(const double step_size);

  void augmentDualResidual(const Robot& robot, const double dtau,
                           Eigen::VectorXd& Cu);

  double residualL1Nrom(const Robot& robot, const double dtau,
                        const Eigen::VectorXd& u);

  double residualSquaredNrom(const Robot& robot, const double dtau, 
                             const Eigen::VectorXd& u);

private:
  unsigned int dimq_, dimv_, dimc_;
  double barrier_;
  Eigen::VectorXd umin_, slack_, dual_, residual_, dslack_, ddual_;
};

} // namespace pdipm
} // namespace idocp


#endif // IDOCP_JOINT_TORQUE_LOWER_LIMITS_PDIPM_HPP_