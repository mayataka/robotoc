#include "constraints/joint_space_constraints.hpp"

#include <assert.h>


namespace idocp {

JointSpaceConstraints::JointSpaceConstraints(const Robot& robot)
  : barrier_(1.0e-04),
    position_upper_limits_(robot),
    position_lower_limits_(robot),
    velocity_upper_limits_(robot),
    velocity_lower_limits_(robot),
    torque_upper_limits_(robot),
    torque_lower_limits_(robot) {
}


void JointSpaceConstraints::setSlackAndDual(const Robot& robot,
                                            const Eigen::VectorXd& q, 
                                            const Eigen::VectorXd& v, 
                                            const Eigen::VectorXd& u) {
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(u.size() == robot.dimv());
  position_upper_limits_.setSlackAndDual(robot, barrier_, q);
  position_lower_limits_.setSlackAndDual(robot, barrier_, q);
  velocity_upper_limits_.setSlackAndDual(robot, barrier_, v);
  velocity_lower_limits_.setSlackAndDual(robot, barrier_, v);
  torque_upper_limits_.setSlackAndDual(robot, barrier_, u);
  torque_lower_limits_.setSlackAndDual(robot, barrier_, u);
}


void JointSpaceConstraints::linearizeConstraint(const Robot& robot, 
                                                const Eigen::VectorXd& q, 
                                                const Eigen::VectorXd& v, 
                                                const Eigen::VectorXd& u) {
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(u.size() == robot.dimv());
  position_upper_limits_.linearizeConstraint(robot, barrier_, q);
  position_lower_limits_.linearizeConstraint(robot, barrier_, q);
  velocity_upper_limits_.linearizeConstraint(robot, barrier_, v);
  velocity_lower_limits_.linearizeConstraint(robot, barrier_, v);
  torque_upper_limits_.linearizeConstraint(robot, barrier_, u);
  torque_lower_limits_.linearizeConstraint(robot, barrier_, u);
}


void JointSpaceConstraints::updateSlackAndDual(const Robot& robot, 
                                               const Eigen::MatrixXd& du_dq,
                                               const Eigen::MatrixXd& du_dv,
                                               const Eigen::MatrixXd& du_da,
                                               const Eigen::VectorXd& dq,
                                               const Eigen::VectorXd& dv, 
                                               const Eigen::VectorXd& da) {
  assert(du_dq.rows() == robot.dimv());
  assert(du_dq.cols() == robot.dimv());
  assert(du_dv.rows() == robot.dimv());
  assert(dq.size() == robot.dimv());
  assert(dv.size() == robot.dimv());
  assert(du.size() == robot.dimv());
  position_upper_limits_.updateSlackAndDual(robot, dq);
  position_lower_limits_.updateSlackAndDual(robot, dq);
  velocity_upper_limits_.updateSlackAndDual(robot, dv);
  velocity_lower_limits_.updateSlackAndDual(robot, dv);
  torque_upper_limits_.updateSlackAndDual(robot, du_dq, du_dv, du_da,
                                          dq, dv, da);
  torque_lower_limits_.updateSlackAndDual(robot, du_dq, du_dv, du_da,
                                          dq, dv, da);
}


void JointSpaceConstraints::condenseSlackAndDual(const Robot& robot, 
                                                 const Eigen::MatrixXd& du_dq, 
                                                 const Eigen::MatrixXd& du_dv,
                                                 const Eigen::MatrixXd& du_da, 
                                                 Eigen::MatrixXd& Cqq, 
                                                 Eigen::MatrixXd& Cqv, 
                                                 Eigen::MatrixXd& Cqa, 
                                                 Eigen::MatrixXd& Cvv, 
                                                 Eigen::MatrixXd& Cva, 
                                                 Eigen::MatrixXd& Caa, 
                                                 Eigen::VectorXd& Cq, 
                                                 Eigen::VectorXd& Cv, 
                                                 Eigen::VectorXd& Ca) {
  assert(du_dq.rows() == robot.dimv());
  assert(du_dq.cols() == robot.dimv());
  assert(du_dv.rows() == robot.dimv());
  assert(Cqq.rows() == robot.dimv());
  assert(Cqq.cols() == robot.dimv());
  assert(Cqv.rows() == robot.dimv());
  assert(Cqv.cols() == robot.dimv());
  assert(Cqa.rows() == robot.dimv());
  assert(Cqa.cols() == robot.dimv());
  assert(Cvv.rows() == robot.dimv());
  assert(Cvv.cols() == robot.dimv());
  assert(Cva.rows() == robot.dimv());
  assert(Cva.cols() == robot.dimv());
  assert(Caa.rows() == robot.dimv());
  assert(Caa.cols() == robot.dimv());
  assert(Cq.size() == robot.dimv());
  assert(Cv.size() == robot.dimv());
  assert(Ca.size() == robot.dimv());
  position_upper_limits_.condenseSlackAndDual(robot, Cqq, Cq);
  position_lower_limits_.condenseSlackAndDual(robot, Cqq, Cq);
  velocity_upper_limits_.condenseSlackAndDual(robot, Cvv, Cv);
  velocity_lower_limits_.condenseSlackAndDual(robot, Cvv, Cv);
  torque_upper_limits_.condenseSlackAndDual(robot, du_dq, du_dv, 
                                            du_da, Cqq, Cqv, Cqa, Cvv, Cva, Caa,
                                            Cq, Cv, Ca);
  torque_lower_limits_.condenseSlackAndDual(robot, du_dq, du_dv, 
                                            du_da, Cqq, Cqv, Cqa, Cvv, Cva, Caa,
                                            Cq, Cv, Ca);
}


void JointSpaceConstraints::augmentDualResidual(const Robot& robot, 
                                                const Eigen::MatrixXd& du_dq, 
                                                const Eigen::MatrixXd& du_dv, 
                                                const Eigen::MatrixXd& du_da, 
                                                Eigen::VectorXd& Cq, 
                                                Eigen::VectorXd& Cv, 
                                                Eigen::VectorXd& Ca) {
  assert(du_dq.rows() == robot.dimv());
  assert(du_dq.cols() == robot.dimv());
  assert(du_dv.rows() == robot.dimv());
  assert(du_dv.cols() == robot.dimv());
  assert(du_da.rows() == robot.dimv());
  assert(du_da.cols() == robot.dimv());
  assert(Cq.size() == robot.dimv());
  assert(Cv.size() == robot.dimv());
  assert(Ca.size() == robot.dimv());
  position_upper_limits_.augmentDualResidual(robot, Cq);
  position_lower_limits_.augmentDualResidual(robot, Cq);
  velocity_upper_limits_.augmentDualResidual(robot, Cv);
  velocity_lower_limits_.augmentDualResidual(robot, Cv);
  torque_upper_limits_.augmentDualResidual(robot, du_dq, du_dv, du_da,  
                                           Cq, Cv, Ca);
  torque_lower_limits_.augmentDualResidual(robot, du_dq, du_dv, du_da,  
                                           Cq, Cv, Ca);
}


double JointSpaceConstraints::squaredConstraintsResidualNrom(
    const Robot& robot, const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
    const Eigen::VectorXd& u) {
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(u.size() == robot.dimv());
  double error = 0;
  error += position_upper_limits_.squaredConstraintsResidualNrom(robot, q);
  error += position_lower_limits_.squaredConstraintsResidualNrom(robot, q);
  error += velocity_upper_limits_.squaredConstraintsResidualNrom(robot, v);
  error += velocity_lower_limits_.squaredConstraintsResidualNrom(robot, v);
  error += torque_upper_limits_.squaredConstraintsResidualNrom(robot, u);
  error += torque_lower_limits_.squaredConstraintsResidualNrom(robot, u);
  return error;
}

} // namespace idocp