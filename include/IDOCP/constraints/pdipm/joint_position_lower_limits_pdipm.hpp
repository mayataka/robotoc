#ifndef IDOCP_CONSTRAINTS_PDIPM_JOINT_POSITION_LOWER_LIMITS_HPP_
#define IDOCP_CONSTRAINTS_PDIPM_JOINT_POSITION_LOWER_LIMITS_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {
namespace pdipm {

class JointPositionLowerLimits {
public:
  JointPositionLowerLimits(const Robot& robot, const double barrier);

  // Use default copy constructor.
  JointPositionLowerLimits(const JointPositionLowerLimits&) = default;

  // Use default copy operator.
  JointPositionLowerLimits& operator=(const JointPositionLowerLimits&) = default;

  bool isFeasible(const Robot& robot, const Eigen::VectorXd& q);

  void setSlackAndDual(const Robot& robot, const double dtau, 
                       const Eigen::VectorXd& q);

  void condenseSlackAndDual(const Robot& robot, const double dtau, 
                            const Eigen::VectorXd& q, Eigen::MatrixXd& Cqq, 
                            Eigen::VectorXd& Cq);

  void computeSlackAndDualDirection(const Robot& robot, const double dtau,
                                    const Eigen::VectorXd& dq);

  double maxSlackStepSize(const double margin_rate);

  double maxDualStepSize(const double margin_rate);

  void updateSlack(const double step_size);

  void updateDual(const double step_size);

  double costSlackBarrier();

  double costSlackBarrier(const double step_size);

  void augmentDualResidual(const Robot& robot, const double dtau, 
                           Eigen::VectorXd& Cq);

  double residualL1Nrom(const Robot& robot, const double dtau,
                        const Eigen::VectorXd& q);

  double residualSquaredNrom(const Robot& robot, const double dtau, 
                             const Eigen::VectorXd& q);

private:
  int dimq_, dimv_, dimc_;
  double barrier_;
  Eigen::VectorXd qmin_, slack_, dual_, residual_, duality_, dslack_, ddual_;
};

} // namespace pdipm
} // namespace idocp


#endif // IDOCP_CONSTRAINTS_PDIPM_JOINT_POSITION_LOWER_LIMITS_HPP_