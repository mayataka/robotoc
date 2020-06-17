#ifndef IDOCP_JOINT_POSITION_LOWER_LIMITS_PDIPM_HPP_
#define IDOCP_JOINT_POSITION_LOWER_LIMITS_PDIPM_HPP_

#include <utility>

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

  std::pair<double, double> computeDirectionAndMaxStepSize(
      const Robot& robot, const double fraction_to_boundary_rate, 
      const double dtau, const Eigen::VectorXd& dq);

  void updateSlack(const double step_size);

  void updateDual(const double step_size);

  double slackBarrier();

  double slackBarrier(const double step_size);

  void augmentDualResidual(const Robot& robot, const double dtau, 
                           Eigen::VectorXd& Cq);

  double residualL1Nrom(const Robot& robot, const double dtau,
                        const Eigen::VectorXd& q);

  double residualSquaredNrom(const Robot& robot, const double dtau, 
                             const Eigen::VectorXd& q);

private:
  unsigned int dimq_, dimv_, dimc_;
  double barrier_;
  Eigen::VectorXd qmin_, slack_, dual_, residual_, dslack_, ddual_;
};

} // namespace pdipm
} // namespace idocp


#endif // IDOCP_JOINT_POSITION_LOWER_LIMITS_PDIPM_HPP_