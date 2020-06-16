#include "constraints/pdipm/joint_space_constraints_pdipm.hpp"

#include <assert.h>


namespace idocp {
namespace pdipm {

JointSpaceConstraints::JointSpaceConstraints(const Robot& robot)
  : barrier_(1.0e-04),
    step_size_reduction_rate_(0.995),
    position_upper_limits_(robot, barrier_),
    position_lower_limits_(robot, barrier_),
    velocity_upper_limits_(robot, barrier_),
    velocity_lower_limits_(robot, barrier_),
    torque_upper_limits_(robot, barrier_),
    torque_lower_limits_(robot, barrier_) {
  assert(barrier_ > 0);
  assert(step_size_reduction_rate_ <= 1);
  assert(step_size_reduction_rate_ > 0);
}


void JointSpaceConstraints::setSlackAndDual(const Robot& robot, 
                                            const double dtau, 
                                            const Eigen::VectorXd& q, 
                                            const Eigen::VectorXd& v, 
                                            const Eigen::VectorXd& a, 
                                            const Eigen::VectorXd& u) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(a.size() == robot.dimv());
  assert(u.size() == robot.dimv());
  position_upper_limits_.setSlackAndDual(robot, dtau, q);
  position_lower_limits_.setSlackAndDual(robot, dtau, q);
  velocity_upper_limits_.setSlackAndDual(robot, dtau, v);
  velocity_lower_limits_.setSlackAndDual(robot, dtau, v);
  torque_upper_limits_.setSlackAndDual(robot, dtau, u);
}



void JointSpaceConstraints::condenseSlackAndDual(const Robot& robot, 
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
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(a.size() == robot.dimv());
  position_upper_limits_.condenseSlackAndDual(robot, dtau, q, Cqq, Cq);
  position_lower_limits_.condenseSlackAndDual(robot, dtau, q, Cqq, Cq);
  velocity_upper_limits_.condenseSlackAndDual(robot, dtau, v, Cvv, Cv);
  velocity_lower_limits_.condenseSlackAndDual(robot, dtau, v, Cvv, Cv);
}

void JointSpaceConstraints::condenseSlackAndDual(const Robot& robot, 
                                                 const double dtau, 
                                                 const Eigen::VectorXd& u, 
                                                 Eigen::MatrixXd& Cuu, 
                                                 Eigen::VectorXd& Cu) {
  assert(dtau > 0);
  assert(u.size() == robot.dimq());
  torque_upper_limits_.condenseSlackAndDual(robot, dtau, u, Cuu, Cu);
  torque_lower_limits_.condenseSlackAndDual(robot, dtau, u, Cuu, Cu);
}


double JointSpaceConstraints::computeDirectionAndMaxStepSize(
    const Robot& robot, const double dtau, const Eigen::MatrixXd& dq, 
    const Eigen::VectorXd& dv, const Eigen::MatrixXd& da, 
    const Eigen::VectorXd& du) {
  assert(dtau > 0);
  assert(dq.size() == robot.dimv());
  assert(dv.size() == robot.dimv());
  assert(da.size() == robot.dimv());
  assert(du.size() == robot.dimv());
  const double position_upper_step 
      = position_upper_limits_.computeDirectionAndMaxStepSize(robot, dtau, dq);
  const double position_lower_step 
      = position_lower_limits_.computeDirectionAndMaxStepSize(robot, dtau, dq);
  const double velocity_upper_step 
      = velocity_upper_limits_.computeDirectionAndMaxStepSize(robot, dtau, dv);
  const double velocity_lower_step 
      = velocity_lower_limits_.computeDirectionAndMaxStepSize(robot, dtau, dv);
  const double torque_upper_step 
      = torque_upper_limits_.computeDirectionAndMaxStepSize(robot, dtau, du);
  const double torque_lower_step 
      = torque_lower_limits_.computeDirectionAndMaxStepSize(robot, dtau, du);
  const double max_step_size =  std::min({position_upper_step, 
                                          position_lower_step, 
                                          velocity_upper_step, 
                                          velocity_lower_step, 
                                          torque_upper_step, 
                                          torque_lower_step});
  if (max_step_size < 1) {
    return max_step_size * step_size_reduction_rate_;
  }
  else {
    return 1;
  }
}


void JointSpaceConstraints::updateSlackAndDual(const Robot& robot, 
                                               const double step_size) {
  assert(step_size > 0);
  position_upper_limits_.updateSlackAndDual(robot, step_size);
  position_lower_limits_.updateSlackAndDual(robot, step_size);
  velocity_upper_limits_.updateSlackAndDual(robot, step_size);
  velocity_lower_limits_.updateSlackAndDual(robot, step_size);
  torque_upper_limits_.updateSlackAndDual(robot, step_size);
  torque_lower_limits_.updateSlackAndDual(robot, step_size);
}


void JointSpaceConstraints::augmentDualResidual(const Robot& robot, 
                                                const double dtau, 
                                                Eigen::VectorXd& Cq,
                                                Eigen::VectorXd& Cv, 
                                                Eigen::VectorXd& Ca, 
                                                Eigen::VectorXd& Cu) {
  assert(dtau > 0);
  assert(Cq.size() == robot.dimv());
  assert(Cv.size() == robot.dimv());
  assert(Ca.size() == robot.dimv());
  assert(Cu.size() == robot.dimv());
  position_upper_limits_.augmentDualResidual(robot, dtau, Cq);
  position_lower_limits_.augmentDualResidual(robot, dtau, Cq);
  velocity_upper_limits_.augmentDualResidual(robot, dtau, Cv);
  velocity_lower_limits_.augmentDualResidual(robot, dtau, Cv);
  torque_upper_limits_.augmentDualResidual(robot, dtau, Cu);
  torque_lower_limits_.augmentDualResidual(robot, dtau, Cu);
}


double JointSpaceConstraints::squaredConstraintsResidualNrom(
    const Robot& robot, const double dtau, const Eigen::VectorXd& q, 
    const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
    const Eigen::VectorXd& u) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(a.size() == robot.dimv());
  assert(u.size() == robot.dimv());
  double norm = 0;
  norm += position_upper_limits_.squaredConstraintsResidualNrom(robot, dtau, q);
  norm += position_lower_limits_.squaredConstraintsResidualNrom(robot, dtau, q);
  norm += velocity_upper_limits_.squaredConstraintsResidualNrom(robot, dtau, v);
  norm += velocity_lower_limits_.squaredConstraintsResidualNrom(robot, dtau, v);
  norm += torque_upper_limits_.squaredConstraintsResidualNrom(robot, dtau, u);
  norm += torque_lower_limits_.squaredConstraintsResidualNrom(robot, dtau, u);
  return norm;
}

} // namespace pdipm
} // namespace idocp