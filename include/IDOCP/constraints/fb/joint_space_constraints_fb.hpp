#ifndef IDOCP_CONSTRAINTS_FB_JOINT_SPACE_CONSTRAINTS_HPP_
#define IDOCP_CONSTRAINTS_FB_JOINT_SPACE_CONSTRAINTS_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"
#include "constraints/fb/joint_position_upper_limits_fb.hpp"
#include "constraints/fb/joint_position_lower_limits_fb.hpp"
#include "constraints/fb/joint_velocity_upper_limits_fb.hpp"
#include "constraints/fb/joint_velocity_lower_limits_fb.hpp"
#include "constraints/fb/joint_torque_upper_limits_fb.hpp"
#include "constraints/fb/joint_torque_lower_limits_fb.hpp"


namespace idocp {
namespace fb {

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

  void augmentDualResidual(const Robot& robot, const double dtau, 
                           Eigen::VectorXd& Cu);

  void augmentDualResidual(const Robot& robot, const double dtau, 
                           Eigen::VectorXd& Cq, Eigen::VectorXd& Cv, 
                           Eigen::VectorXd& Ca);

  void condenseSlackAndDual(const Robot& robot, const double dtau, 
                            const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                            const Eigen::VectorXd& a, Eigen::MatrixXd& Cqq, 
                            Eigen::MatrixXd& Cvv, Eigen::MatrixXd& Caa,  
                            Eigen::VectorXd& Cq, Eigen::VectorXd& Cv, 
                            Eigen::VectorXd& Ca);

  void condenseSlackAndDual(const Robot& robot, const double dtau, 
                            const Eigen::VectorXd& u, Eigen::MatrixXd& Cuu, 
                            Eigen::VectorXd& Cu);

  void computeSlackAndDualDirection(const Robot& robot, const double dtau,
                                    const Eigen::VectorXd& dq,
                                    const Eigen::VectorXd& dv,
                                    const Eigen::VectorXd& da,
                                    const Eigen::VectorXd& du);

  double maxSlackStepSize();

  double maxDualStepSize();

  void updateSlack(const double step_size);

  void updateDual(const double step_size);

  double costSlackBarrier();

  double costSlackBarrier(const double step_size);

  double residualL1Nrom(const Robot& robot, const double dtau,
                        const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                        const Eigen::VectorXd& a, const Eigen::VectorXd& u);

  double residualSquaredNrom(const Robot& robot, const double dtau,
                             const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                             const Eigen::VectorXd& a, const Eigen::VectorXd& u);

private:
  double barrier_, fraction_to_boundary_margin_;
  fb::JointPositionUpperLimits position_upper_limits_;
  fb::JointPositionLowerLimits position_lower_limits_;
  fb::JointVelocityUpperLimits velocity_upper_limits_;
  fb::JointVelocityLowerLimits velocity_lower_limits_;
  fb::JointTorqueUpperLimits torque_upper_limits_;
  fb::JointTorqueLowerLimits torque_lower_limits_;
};

} // namespace fb
} // namespace idocp


#endif // IDOCP_CONSTRAINTS_FB_JOINT_SPACE_CONSTRAINTS_HPP_