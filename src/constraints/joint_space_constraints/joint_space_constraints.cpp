#include "constraints/joint_space_constraints/joint_space_constraints.hpp"

#include <assert.h>


namespace idocp {
namespace pdipm {

JointSpaceConstraints::JointSpaceConstraints(const Robot& robot)
  : time_step_(2),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    barrier_(1.0e-04),
    fraction_to_boundary_margin_(0.995),
    position_upper_limits_(robot, robot.upperJointPositionLimit(), barrier_),
    velocity_upper_limits_(robot, robot.jointVelocityLimit(), barrier_),
    torque_upper_limits_(robot, robot.jointEffortLimit(), barrier_),
    position_lower_limits_(robot, robot.lowerJointPositionLimit(), barrier_),
    velocity_lower_limits_(robot, -robot.jointVelocityLimit(), barrier_),
    torque_lower_limits_(robot, -robot.jointEffortLimit(), barrier_) {
  assert(barrier_ > 0);
  assert(fraction_to_boundary_margin_ <= 1);
  assert(fraction_to_boundary_margin_ > 0);
}


void JointSpaceConstraints::setTimeStep(const int time_step) {
  assert(time_step >= 0);
  time_step_ = time_step;
}


bool JointSpaceConstraints::isFeasible(const Eigen::VectorXd& q, 
                                       const Eigen::VectorXd& v, 
                                       const Eigen::VectorXd& a, 
                                       const Eigen::VectorXd& u) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(u.size() == dimv_);
  if (time_step_ >= 2) {
    if (!position_upper_limits_.isFeasible(q)) {
      return false;
    }
    if (!position_lower_limits_.isFeasible(q)) {
      return false;
    }
  }
  if (time_step_ >= 1) {
    if (!velocity_upper_limits_.isFeasible(v)) {
      return false;
    }
    if (!velocity_lower_limits_.isFeasible(v)) {
      return false;
    }
  }
  if (!torque_upper_limits_.isFeasible(u)) {
    return false;
  }
  if (!torque_lower_limits_.isFeasible(u)) {
    return false;
  }
  return true;
}


void JointSpaceConstraints::setSlackAndDual(const double dtau, 
                                            const Eigen::VectorXd& q, 
                                            const Eigen::VectorXd& v, 
                                            const Eigen::VectorXd& a, 
                                            const Eigen::VectorXd& u) {
  assert(dtau > 0);
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(u.size() == dimv_);
  if (time_step_ >= 2) {
    position_upper_limits_.setSlackAndDual(dtau, q);
    position_lower_limits_.setSlackAndDual(dtau, q);
  }
  if (time_step_ >= 1) {
    velocity_upper_limits_.setSlackAndDual(dtau, v);
    velocity_lower_limits_.setSlackAndDual(dtau, v);
  }
  torque_upper_limits_.setSlackAndDual(dtau, u);
  torque_lower_limits_.setSlackAndDual(dtau, u);
}


void JointSpaceConstraints::augmentDualResidual(const double dtau, 
                                                Eigen::VectorXd& Cu) {
  assert(dtau > 0);
  assert(Cu.size() == dimv_);
  torque_upper_limits_.augmentDualResidual(dtau, Cu);
  torque_lower_limits_.augmentDualResidual(dtau, Cu);
}


void JointSpaceConstraints::augmentDualResidual(const double dtau, 
                                                Eigen::VectorXd& Cq,
                                                Eigen::VectorXd& Cv, 
                                                Eigen::VectorXd& Ca) {
  assert(dtau > 0);
  assert(Cq.size() == dimv_);
  assert(Cv.size() == dimv_);
  assert(Ca.size() == dimv_);
  if (time_step_ >= 2) {
    position_upper_limits_.augmentDualResidual(dtau, Cq);
    position_lower_limits_.augmentDualResidual(dtau, Cq);
  }
  if (time_step_ >= 1) {
    velocity_upper_limits_.augmentDualResidual(dtau, Cv);
    velocity_lower_limits_.augmentDualResidual(dtau, Cv);
  }
}


void JointSpaceConstraints::condenseSlackAndDual(const double dtau, 
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
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(Cqq.rows() == dimv_);
  assert(Cqq.cols() == dimv_);
  assert(Cvv.rows() == dimv_);
  assert(Cvv.cols() == dimv_);
  assert(Caa.rows() == dimv_);
  assert(Caa.cols() == dimv_);
  assert(Cq.size() == dimv_);
  assert(Cv.size() == dimv_);
  assert(Ca.size() == dimv_);
  if (time_step_ >= 2) {
    position_upper_limits_.condenseSlackAndDual(dtau, q, Cqq, Cq);
    position_lower_limits_.condenseSlackAndDual(dtau, q, Cqq, Cq);
  }
  if (time_step_ >= 1) {
    velocity_upper_limits_.condenseSlackAndDual(dtau, v, Cvv, Cv);
    velocity_lower_limits_.condenseSlackAndDual(dtau, v, Cvv, Cv);
  }
}


void JointSpaceConstraints::condenseSlackAndDual(const double dtau, 
                                                 const Eigen::VectorXd& u, 
                                                 Eigen::MatrixXd& Cuu, 
                                                 Eigen::VectorXd& Cu) {
  assert(dtau > 0);
  assert(u.size() == dimv_);
  torque_upper_limits_.condenseSlackAndDual(dtau, u, Cuu, Cu);
  torque_lower_limits_.condenseSlackAndDual(dtau, u, Cuu, Cu);
}


void JointSpaceConstraints::computeSlackAndDualDirection(
    const double dtau, const Eigen::VectorXd& dq, const Eigen::VectorXd& dv, 
    const Eigen::VectorXd& da, const Eigen::VectorXd& du) {
  assert(dtau > 0);
  assert(dq.size() == dimv_);
  assert(dv.size() == dimv_);
  assert(da.size() == dimv_);
  assert(du.size() == dimv_);
  if (time_step_ >= 2) {
    position_upper_limits_.computeSlackAndDualDirection(dtau, dq);
    position_lower_limits_.computeSlackAndDualDirection(dtau, dq);
  }
  if (time_step_ >= 1) {
    velocity_upper_limits_.computeSlackAndDualDirection(dtau, dv);
    velocity_lower_limits_.computeSlackAndDualDirection(dtau, dv);
  }
  torque_upper_limits_.computeSlackAndDualDirection(dtau, du);
  torque_lower_limits_.computeSlackAndDualDirection(dtau, du);
}


double JointSpaceConstraints::maxSlackStepSize() {
  assert(fraction_to_boundary_margin_ > 0);
  assert(fraction_to_boundary_margin_ <= 1);
  if (time_step_ >= 2) {
    const double size_position_upper_limit 
        = position_upper_limits_.maxSlackStepSize(fraction_to_boundary_margin_);
    const double size_position_lower_limit 
        = position_lower_limits_.maxSlackStepSize(fraction_to_boundary_margin_);
    const double size_velocity_upper_limit 
        = velocity_upper_limits_.maxSlackStepSize(fraction_to_boundary_margin_);
    const double size_velocity_lower_limit 
        = velocity_lower_limits_.maxSlackStepSize(fraction_to_boundary_margin_);
    const double size_torque_upper_limit 
        = torque_upper_limits_.maxSlackStepSize(fraction_to_boundary_margin_);
    const double size_torque_lower_limit 
        = torque_lower_limits_.maxSlackStepSize(fraction_to_boundary_margin_);
    return std::min({size_position_upper_limit, size_position_lower_limit, 
                     size_velocity_upper_limit, size_velocity_lower_limit,
                     size_torque_upper_limit, size_torque_lower_limit});
  }
  else if (time_step_ >= 1) {
    const double size_velocity_upper_limit 
        = velocity_upper_limits_.maxSlackStepSize(fraction_to_boundary_margin_);
    const double size_velocity_lower_limit 
        = velocity_lower_limits_.maxSlackStepSize(fraction_to_boundary_margin_);
    const double size_torque_upper_limit 
        = torque_upper_limits_.maxSlackStepSize(fraction_to_boundary_margin_);
    const double size_torque_lower_limit 
        = torque_lower_limits_.maxSlackStepSize(fraction_to_boundary_margin_);
    return std::min({size_velocity_upper_limit, size_velocity_lower_limit,
                     size_torque_upper_limit, size_torque_lower_limit});
  }
  const double size_torque_upper_limit 
      = torque_upper_limits_.maxSlackStepSize(fraction_to_boundary_margin_);
  const double size_torque_lower_limit 
      = torque_lower_limits_.maxSlackStepSize(fraction_to_boundary_margin_);
  return std::min({size_torque_upper_limit, size_torque_lower_limit});
}


double JointSpaceConstraints::maxDualStepSize() {
  assert(fraction_to_boundary_margin_ > 0);
  assert(fraction_to_boundary_margin_ <= 1);
  if (time_step_ >= 2) {
    const double size_position_upper_limit 
        = position_upper_limits_.maxDualStepSize(fraction_to_boundary_margin_);
    const double size_position_lower_limit 
        = position_lower_limits_.maxDualStepSize(fraction_to_boundary_margin_);
    const double size_velocity_upper_limit 
        = velocity_upper_limits_.maxDualStepSize(fraction_to_boundary_margin_);
    const double size_velocity_lower_limit 
        = velocity_lower_limits_.maxDualStepSize(fraction_to_boundary_margin_);
    const double size_torque_upper_limit 
        = torque_upper_limits_.maxDualStepSize(fraction_to_boundary_margin_);
    const double size_torque_lower_limit 
        = torque_lower_limits_.maxDualStepSize(fraction_to_boundary_margin_);
    return std::min({size_position_upper_limit, size_position_lower_limit, 
                     size_velocity_upper_limit, size_velocity_lower_limit,
                     size_torque_upper_limit, size_torque_lower_limit});
  }
  else if (time_step_ >= 1) {
    const double size_velocity_upper_limit 
        = velocity_upper_limits_.maxDualStepSize(fraction_to_boundary_margin_);
    const double size_velocity_lower_limit 
        = velocity_lower_limits_.maxDualStepSize(fraction_to_boundary_margin_);
    const double size_torque_upper_limit 
        = torque_upper_limits_.maxDualStepSize(fraction_to_boundary_margin_);
    const double size_torque_lower_limit 
        = torque_lower_limits_.maxDualStepSize(fraction_to_boundary_margin_);
    return std::min({size_velocity_upper_limit, size_velocity_lower_limit,
                     size_torque_upper_limit, size_torque_lower_limit});
  }
  const double size_torque_upper_limit 
      = torque_upper_limits_.maxDualStepSize(fraction_to_boundary_margin_);
  const double size_torque_lower_limit 
      = torque_lower_limits_.maxDualStepSize(fraction_to_boundary_margin_);
  return std::min({size_torque_upper_limit, size_torque_lower_limit});
}


void JointSpaceConstraints::updateSlack(const double step_size) {
  assert(step_size > 0);
  if (time_step_ >= 2) {
    position_upper_limits_.updateSlack(step_size);
    position_lower_limits_.updateSlack(step_size);
  }
  if (time_step_ >= 1) {
    velocity_upper_limits_.updateSlack(step_size);
    velocity_lower_limits_.updateSlack(step_size);
  }
  torque_upper_limits_.updateSlack(step_size);
  torque_lower_limits_.updateSlack(step_size);
}


void JointSpaceConstraints::updateDual(const double step_size) {
  assert(step_size > 0);
  if (time_step_ >= 2) {
    position_upper_limits_.updateDual(step_size);
    position_lower_limits_.updateDual(step_size);
  }
  if (time_step_ >= 1) {
    velocity_upper_limits_.updateDual(step_size);
    velocity_lower_limits_.updateDual(step_size);
  }
  torque_upper_limits_.updateDual(step_size);
  torque_lower_limits_.updateDual(step_size);
}


double JointSpaceConstraints::costSlackBarrier() {
  double cost = 0;
  if (time_step_ >= 2) {
    cost += position_upper_limits_.costSlackBarrier();
    cost += position_lower_limits_.costSlackBarrier();
  }
  if (time_step_ >= 1) {
    cost += velocity_upper_limits_.costSlackBarrier();
    cost += velocity_lower_limits_.costSlackBarrier();
  }
  cost += torque_upper_limits_.costSlackBarrier();
  cost += torque_lower_limits_.costSlackBarrier();
  return barrier_ * cost;
}


double JointSpaceConstraints::costSlackBarrier(const double step_size) {
  double cost = 0;
  if (time_step_ >= 2) {
    cost += position_upper_limits_.costSlackBarrier(step_size);
    cost += position_lower_limits_.costSlackBarrier(step_size);
  }
  if (time_step_ >= 1) {
    cost += velocity_upper_limits_.costSlackBarrier(step_size);
    cost += velocity_lower_limits_.costSlackBarrier(step_size);
  }
  cost += torque_upper_limits_.costSlackBarrier(step_size);
  cost += torque_lower_limits_.costSlackBarrier(step_size);
  return barrier_ * cost;
}


double JointSpaceConstraints::residualL1Nrom(const double dtau, 
                                             const Eigen::VectorXd& q, 
                                             const Eigen::VectorXd& v, 
                                             const Eigen::VectorXd& a, 
                                             const Eigen::VectorXd& u) {
  assert(dtau > 0);
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(u.size() == dimv_);
  double norm = 0;
  if (time_step_ >= 2) {
    norm += position_upper_limits_.residualL1Nrom(dtau, q);
    norm += position_lower_limits_.residualL1Nrom(dtau, q);
  }
  if (time_step_ >= 1) {
    norm += velocity_upper_limits_.residualL1Nrom(dtau, v);
    norm += velocity_lower_limits_.residualL1Nrom(dtau, v);
  }
  norm += torque_upper_limits_.residualL1Nrom(dtau, u);
  norm += torque_lower_limits_.residualL1Nrom(dtau, u);
  return norm;
}


double JointSpaceConstraints::residualSquaredNrom(const double dtau, 
                                                  const Eigen::VectorXd& q, 
                                                  const Eigen::VectorXd& v, 
                                                  const Eigen::VectorXd& a, 
                                                  const Eigen::VectorXd& u) {
  assert(dtau > 0);
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(u.size() == dimv_);
  double norm = 0;
  if (time_step_ >= 2) {
    norm += position_upper_limits_.residualSquaredNrom(dtau, q);
    norm += position_lower_limits_.residualSquaredNrom(dtau, q);
  }
  if (time_step_ >= 1) {
    norm += velocity_upper_limits_.residualSquaredNrom(dtau, v);
    norm += velocity_lower_limits_.residualSquaredNrom(dtau, v);
  }
  norm += torque_upper_limits_.residualSquaredNrom(dtau, u);
  norm += torque_lower_limits_.residualSquaredNrom(dtau, u);
  return norm;
}

} // namespace pdipm
} // namespace idocp