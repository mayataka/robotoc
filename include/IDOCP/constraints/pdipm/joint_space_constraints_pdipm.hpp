#ifndef IDOCP_JOINT_SPACE_CONSTRAINTS_PDIPM_HPP_
#define IDOCP_JOINT_SPACE_CONSTRAINTS_PDIPM_HPP_

#include <utility> 

#include "Eigen/Core"

#include "robot/robot.hpp"
#include "constraints/pdipm/joint_position_upper_limits_pdipm.hpp"
#include "constraints/pdipm/joint_position_lower_limits_pdipm.hpp"
#include "constraints/pdipm/joint_velocity_upper_limits_pdipm.hpp"
#include "constraints/pdipm/joint_velocity_lower_limits_pdipm.hpp"
#include "constraints/pdipm/joint_torque_upper_limits_pdipm.hpp"
#include "constraints/pdipm/joint_torque_lower_limits_pdipm.hpp"


namespace idocp {
namespace pdipm {

class JointSpaceConstraints {
public:
  JointSpaceConstraints(const Robot& robot);

  // Use default copy constructor.
  JointSpaceConstraints(const JointSpaceConstraints&) = default;

  // Use default copy operator.
  JointSpaceConstraints& operator=(const JointSpaceConstraints&) = default;

  bool isFeasible(const Robot& robot, const Eigen::VectorXd& q, 
                  const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                  const Eigen::VectorXd& u);

  void setSlackAndDual(const Robot& robot, const double dtau, 
                       const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                       const Eigen::VectorXd& a, const Eigen::VectorXd& u);

  void condenseSlackAndDual(const Robot& robot, const double dtau, 
                            const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                            const Eigen::VectorXd& a, Eigen::MatrixXd& Cqq, 
                            Eigen::MatrixXd& Cvv, Eigen::MatrixXd& Caa,  
                            Eigen::VectorXd& Cq, Eigen::VectorXd& Cv, 
                            Eigen::VectorXd& Ca);

  void condenseSlackAndDual(const Robot& robot, const double dtau, 
                            const Eigen::VectorXd& u, Eigen::MatrixXd& Cuu, 
                            Eigen::VectorXd& Cu);

  std::pair<double, double> computeDirectionAndMaxStepSize(
      const Robot& robot, const double dtau, const Eigen::MatrixXd& dq, 
      const Eigen::VectorXd& dv, const Eigen::MatrixXd& da, 
      const Eigen::VectorXd& du);

  void updateSlack(const double step_size);

  void updateDual(const double step_size);

  double slackBarrier();

  double slackBarrier(const double step_size);

  void augmentDualResidual(const Robot& robot, const double dtau, 
                           Eigen::VectorXd& Cq,
                           Eigen::VectorXd& Cv, 
                           Eigen::VectorXd& Ca, 
                           Eigen::VectorXd& Cu);

  double residualL1Nrom(const Robot& robot, const double dtau,
                        const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                        const Eigen::VectorXd& a, const Eigen::VectorXd& u);

  double residualSquaredNrom(const Robot& robot, const double dtau,
                             const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                             const Eigen::VectorXd& a, const Eigen::VectorXd& u);

private:
  double barrier_, fraction_to_boundary_rate_;
  pdipm::JointPositionUpperLimits position_upper_limits_;
  pdipm::JointPositionLowerLimits position_lower_limits_;
  pdipm::JointVelocityUpperLimits velocity_upper_limits_;
  pdipm::JointVelocityLowerLimits velocity_lower_limits_;
  pdipm::JointTorqueUpperLimits torque_upper_limits_;
  pdipm::JointTorqueLowerLimits torque_lower_limits_;
};

} // namespace pdipm
} // namespace idocp


#endif // IDOCP_JOINT_SPACE_CONSTRAINTS_PDIPM_HPP_