#include "constraints/joint_space_soft_constraints.hpp"

#include <assert.h>


namespace idocp {

JointSpaceSoftConstraints::JointSpaceSoftConstraints(const Robot& robot)
  : position_upper_limits_(robot, 1.0e-04),
    position_lower_limits_(robot, 1.0e-04),
    velocity_upper_limits_(robot, 1.0e-04),
    velocity_lower_limits_(robot, 1.0e-04),
    torque_upper_limits_(robot, 1.0e-04),
    torque_lower_limits_(robot, 1.0e-04) {
}


void JointSpaceSoftConstraints::setSlackAndDual(const Robot& robot,
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
  torque_lower_limits_.setSlackAndDual(robot, dtau, u);
}


void JointSpaceSoftConstraints::linearizeConstraint(const Robot& robot, 
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
  position_upper_limits_.linearizeConstraint(robot, dtau, q);
  position_lower_limits_.linearizeConstraint(robot, dtau, q);
  velocity_upper_limits_.linearizeConstraint(robot, dtau, v);
  velocity_lower_limits_.linearizeConstraint(robot, dtau, v);
  torque_upper_limits_.linearizeConstraint(robot, dtau, u);
  torque_lower_limits_.linearizeConstraint(robot, dtau, u);
}


void JointSpaceSoftConstraints::updateSlackAndDual(const Robot& robot, 
                                                   const double dtau,
                                                   const Eigen::VectorXd& dq,
                                                   const Eigen::VectorXd& dv,
                                                   const Eigen::VectorXd& da, 
                                                   const Eigen::VectorXd& du) {
  assert(dtau > 0);
  assert(dq.size() == robot.dimv());
  assert(dv.size() == robot.dimv());
  assert(da.size() == robot.dimv());
  assert(du.size() == robot.dimv());
  position_upper_limits_.updateSlackAndDual(robot, dtau, dq);
  position_lower_limits_.updateSlackAndDual(robot, dtau, dq);
  velocity_upper_limits_.updateSlackAndDual(robot, dtau, dv);
  velocity_lower_limits_.updateSlackAndDual(robot, dtau, dv);
  torque_upper_limits_.updateSlackAndDual(robot, dtau, du);
  torque_lower_limits_.updateSlackAndDual(robot, dtau, du);
}


void JointSpaceSoftConstraints::augmentLuAndLuu(const Robot& robot, 
                                                const double dtau, 
                                                const Eigen::VectorXd& u, 
                                                Eigen::VectorXd& lu, 
                                                Eigen::MatrixXd& luu) {
  assert(dtau > 0);
  assert(u.size() == robot.dimv());
  assert(lu.size() == robot.dimv());
  assert(luu.rows() == robot.dimv());
  assert(luu.cols() == robot.dimv());
  torque_upper_limits_.augmentLuAndLuu(robot, dtau, u, lu, luu);
  torque_lower_limits_.augmentLuAndLuu(robot, dtau, u, lu, luu);
}


void JointSpaceSoftConstraints::condenseSlackAndDual(const Robot& robot, 
                                                     const double dtau,
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
  assert(dtau > 0);
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
  position_upper_limits_.condenseSlackAndDual(robot, dtau, Cqq, Cq);
  position_lower_limits_.condenseSlackAndDual(robot, dtau, Cqq, Cq);
  velocity_upper_limits_.condenseSlackAndDual(robot, dtau, Cvv, Cv);
  velocity_lower_limits_.condenseSlackAndDual(robot, dtau, Cvv, Cv);
  torque_upper_limits_.condenseSlackAndDual(robot, dtau, du_dq, du_dv, 
                                            du_da, Cqq, Cqv, Cqa, Cvv, Cva, Caa,
                                            Cq, Cv, Ca);
  torque_lower_limits_.condenseSlackAndDual(robot, dtau, du_dq, du_dv, 
                                            du_da, Cqq, Cqv, Cqa, Cvv, Cva, Caa,
                                            Cq, Cv, Ca);
}


void JointSpaceSoftConstraints::augmentDualResidual(const Robot& robot, 
                                                    const double dtau,
                                                    const Eigen::MatrixXd& du_dq, 
                                                    const Eigen::MatrixXd& du_dv, 
                                                    const Eigen::MatrixXd& du_da, 
                                                    Eigen::VectorXd& Cq, 
                                                    Eigen::VectorXd& Cv, 
                                                    Eigen::VectorXd& Ca) {
  assert(dtau > 0);
  assert(du_dq.rows() == robot.dimv());
  assert(du_dq.cols() == robot.dimv());
  assert(du_dv.rows() == robot.dimv());
  assert(du_dv.cols() == robot.dimv());
  assert(du_da.rows() == robot.dimv());
  assert(du_da.cols() == robot.dimv());
  assert(Cq.size() == robot.dimv());
  assert(Cv.size() == robot.dimv());
  assert(Ca.size() == robot.dimv());
  position_upper_limits_.augmentDualResidual(robot, dtau, Cq);
  position_lower_limits_.augmentDualResidual(robot, dtau, Cq);
  velocity_upper_limits_.augmentDualResidual(robot, dtau, Cv);
  velocity_lower_limits_.augmentDualResidual(robot, dtau, Cv);
  torque_upper_limits_.augmentDualResidual(robot, dtau, du_dq, du_dv, du_da,  
                                           Cq, Cv, Ca);
  torque_lower_limits_.augmentDualResidual(robot, dtau, du_dq, du_dv, du_da,  
                                           Cq, Cv, Ca);
}


double JointSpaceSoftConstraints::squaredConstraintsResidualNrom(
    const Robot& robot, const double dtau, const Eigen::VectorXd& q, 
    const Eigen::VectorXd& v, const Eigen::VectorXd& u) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(u.size() == robot.dimv());
  double error = 0;
  error += position_upper_limits_.squaredConstraintsResidualNrom(robot, dtau, q);
  error += position_lower_limits_.squaredConstraintsResidualNrom(robot, dtau, q);
  error += velocity_upper_limits_.squaredConstraintsResidualNrom(robot, dtau, v);
  error += velocity_lower_limits_.squaredConstraintsResidualNrom(robot, dtau, v);
  error += torque_upper_limits_.squaredConstraintsResidualNrom(robot, dtau, u);
  error += torque_lower_limits_.squaredConstraintsResidualNrom(robot, dtau, u);
  return error;
}

} // namespace idocp