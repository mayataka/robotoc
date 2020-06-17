#ifndef IDOCP_JOINT_POSITION_UPPER_LIMITS_PDIPM_HPP_
#define IDOCP_JOINT_POSITION_UPPER_LIMITS_PDIPM_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {
namespace pdipm {

class JointPositionUpperLimits {
public:
  JointPositionUpperLimits(const Robot& robot, const double barrier);

  // Use default copy constructor.
  JointPositionUpperLimits(const JointPositionUpperLimits&) = default;

  // Use default copy operator.
  JointPositionUpperLimits& operator=(const JointPositionUpperLimits&) = default;

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
  Eigen::VectorXd qmax_, slack_, dual_, residual_, dslack_, ddual_; 
};

} // namespace pdipm
} // namespace idocp


#endif // IDOCP_JOINT_POSITION_UPPER_LIMITS_PDIPM_HPP_