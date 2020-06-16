#ifndef IDOCP_JOINT_SPACE_CONSTRAINTS_BARRIER_HPP_
#define IDOCP_JOINT_SPACE_CONSTRAINTS_BARRIER_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"
#include "constraints/barrier/joint_position_upper_limits_barrier.hpp"
#include "constraints/barrier/joint_position_lower_limits_barrier.hpp"
#include "constraints/barrier/joint_velocity_upper_limits_barrier.hpp"
#include "constraints/barrier/joint_velocity_lower_limits_barrier.hpp"
#include "constraints/barrier/joint_torque_upper_limits_barrier.hpp"
#include "constraints/barrier/joint_torque_lower_limits_barrier.hpp"


namespace idocp {

class JointSpaceConstraintsBarrier {
public:
  JointSpaceConstraintsBarrier(const Robot& robot);

  // Use default copy constructor.
  JointSpaceConstraintsBarrier(const JointSpaceConstraintsBarrier&) = default;

  // Use default copy operator.
  JointSpaceConstraintsBarrier& operator=(const JointSpaceConstraintsBarrier&) = default;

  void augmentBarrier(const Robot& robot, const double dtau, 
                      const Eigen::VectorXd& u, Eigen::MatrixXd& Cuu, 
                      Eigen::VectorXd& Cu);

  void augmentBarrier(const Robot& robot, const double dtau, 
                      const Eigen::VectorXd& q, const Eigen::VectorXd& v,
                      const Eigen::VectorXd& a, Eigen::MatrixXd& Cqq, 
                      Eigen::MatrixXd& Cvv, Eigen::MatrixXd& Caa, 
                      Eigen::VectorXd& Cq, Eigen::VectorXd& Cv, 
                      Eigen::VectorXd& Ca);

  void augmentBarrierResidual(const Robot& robot, const double dtau, 
                              const Eigen::VectorXd& q, 
                              const Eigen::VectorXd& v,
                              const Eigen::VectorXd& a, 
                              const Eigen::VectorXd& u,
                              Eigen::VectorXd& Cq, Eigen::VectorXd& Cv, 
                              Eigen::VectorXd& Ca, Eigen::VectorXd& Cu);

private:
  double barrier_;
  barrier::JointPositionUpperLimits position_upper_limits_;
  barrier::JointPositionLowerLimits position_lower_limits_;
  barrier::JointVelocityUpperLimits velocity_upper_limits_;
  barrier::JointVelocityLowerLimits velocity_lower_limits_;
  barrier::JointTorqueUpperLimits torque_upper_limits_;
  barrier::JointTorqueLowerLimits torque_lower_limits_;
};

} // namespace idocp


#endif // IDOCP_JOINT_SPACE_CONSTRAINTS_BARRIER_HPP_