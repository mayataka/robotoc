#ifndef IDOCP_CONSTRAINTS_JOINT_TORQUE_UPPER_LIMITS_HPP_
#define IDOCP_CONSTRAINTS_JOINT_TORQUE_UPPER_LIMITS_HPP_

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

  bool isFeasible(const Robot& robot, const Eigen::VectorXd& u);

  void setSlackAndDual(const Robot& robot, const double dtau,
                       const Eigen::VectorXd& u);

  void condenseSlackAndDual(const Robot& robot, const double dtau, 
                            const Eigen::VectorXd& u, Eigen::MatrixXd& Cuu, 
                            Eigen::VectorXd& Cu);

  void computeSlackAndDualDirection(const Robot& robot, const double dtau,
                                    const Eigen::VectorXd& dq);

  double maxSlackStepSize(const double margin_rate);

  double maxDualStepSize(const double margin_rate);

  void updateSlack(const double step_size);

  void updateDual(const double step_size);

  double costSlackBarrier();

  double costSlackBarrier(const double step_size);

  void augmentDualResidual(const Robot& robot, const double dtau,
                           Eigen::VectorXd& Cu);

  double residualL1Nrom(const Robot& robot, const double dtau,
                        const Eigen::VectorXd& u);

  double residualSquaredNrom(const Robot& robot, const double dtau, 
                             const Eigen::VectorXd& u);

private:
  int dimq_, dimv_, dimc_;
  double barrier_;
  Eigen::VectorXd umax_, slack_, dual_, residual_, duality_, dslack_, ddual_;
};

} // namespace pdipm
} // namespace idocp


#endif // IDOCP_CONSTRAINTS_JOINT_TORQUE_UPPER_LIMITS_HPP_