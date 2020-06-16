#include "constraints/barrier/joint_space_constraints_barrier.hpp"

#include <assert.h>


namespace idocp {

JointSpaceConstraintsBarrier::JointSpaceConstraintsBarrier(const Robot& robot)
  : position_upper_limits_(robot, 1.0e-04),
    position_lower_limits_(robot, 1.0e-04),
    velocity_upper_limits_(robot, 1.0e-04),
    velocity_lower_limits_(robot, 1.0e-04),
    torque_upper_limits_(robot, 1.0e-04),
    torque_lower_limits_(robot, 1.0e-04) {
}


void JointSpaceConstraintsBarrier::augmentBarrier(const Robot& robot, 
                                                  const double dtau, 
                                                  const Eigen::VectorXd& u, 
                                                  Eigen::MatrixXd& Cuu, 
                                                  Eigen::VectorXd& Cu) {
  torque_upper_limits_.augmentBarrier(robot, dtau, u, Cuu, Cu);
  torque_lower_limits_.augmentBarrier(robot, dtau, u, Cuu, Cu);
}


void JointSpaceConstraintsBarrier::augmentBarrier(const Robot& robot, 
                                                  const double dtau, 
                                                  const Eigen::VectorXd& q, 
                                                  const Eigen::VectorXd& v, 
                                                  const Eigen::VectorXd& a, 
                                                  Eigen::MatrixXd& Cqq, 
                                                  Eigen::MatrixXd& Cvv, 
                                                  Eigen::MatrixXd& Caa, 
                                                  Eigen::VectorXd& Cq, 
                                                  Eigen::VectorXd& Cv, 
                                                  Eigen::VectorXd& Ca) {
  position_upper_limits_.augmentBarrier(robot, dtau, q, Cqq, Cq);
  position_lower_limits_.augmentBarrier(robot, dtau, q, Cqq, Cq);
  velocity_upper_limits_.augmentBarrier(robot, dtau, v, Cvv, Cv);
  velocity_lower_limits_.augmentBarrier(robot, dtau, v, Cvv, Cv);
}


void JointSpaceConstraintsBarrier::augmentBarrierResidual(
    const Robot& robot, const double dtau, const Eigen::VectorXd& q, 
    const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
    const Eigen::VectorXd& u, Eigen::VectorXd& Cq, Eigen::VectorXd& Cv, 
    Eigen::VectorXd& Ca, Eigen::VectorXd& Cu) {
  position_upper_limits_.augmentBarrierResidual(robot, dtau, q, Cq);
  position_lower_limits_.augmentBarrierResidual(robot, dtau, q, Cq);
  velocity_upper_limits_.augmentBarrierResidual(robot, dtau, v, Cv);
  velocity_lower_limits_.augmentBarrierResidual(robot, dtau, v, Cv);
  torque_upper_limits_.augmentBarrierResidual(robot, dtau, u, Cu);
  torque_lower_limits_.augmentBarrierResidual(robot, dtau, u, Cu);
}

} // namespace idocp