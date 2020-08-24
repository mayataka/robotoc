#ifndef IDOCP_CONSTRAINTS_JOINT_SPACE_CONSTRAINTS_HPP_
#define IDOCP_CONSTRAINTS_JOINT_SPACE_CONSTRAINTS_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "joint_variables_upper_limits.hpp"
#include "joint_variables_lower_limits.hpp"


namespace idocp {
namespace pdipmtest {

class JointSpaceConstraints {
public:
  JointSpaceConstraints(const Robot& robot);

  JointSpaceConstraints();

  ~JointSpaceConstraints();

  // Use default copy constructor.
  JointSpaceConstraints(const JointSpaceConstraints&) = default;

  // Use default copy operator.
  JointSpaceConstraints& operator=(const JointSpaceConstraints&) = default;

  // Use default move constructor.
  JointSpaceConstraints(JointSpaceConstraints&&) noexcept = default;

  // Use default move assign operator.
  JointSpaceConstraints& operator=(JointSpaceConstraints&&) noexcept = default;

  void setTimeStep(const int time_step);

  bool isFeasible(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u);

  void setSlackAndDual(const double dtau, const Eigen::VectorXd& q, 
                       const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                       const Eigen::VectorXd& u);

  void augmentDualResidual(const double dtau, Eigen::VectorXd& Cu);

  void augmentDualResidual(const double dtau, Eigen::VectorXd& Cq, 
                           Eigen::VectorXd& Cv, Eigen::VectorXd& Ca);

  void condenseSlackAndDual(const double dtau, const Eigen::VectorXd& u, 
                            Eigen::MatrixXd& Cuu, Eigen::VectorXd& Cu);

  void condenseSlackAndDual(const double dtau, const Eigen::VectorXd& q, 
                            const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                            Eigen::MatrixXd& Cqq, Eigen::MatrixXd& Cvv, 
                            Eigen::MatrixXd& Caa,  Eigen::VectorXd& Cq, 
                            Eigen::VectorXd& Cv, Eigen::VectorXd& Ca);

  void computeSlackAndDualDirection(const double dtau,
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

  double residualL1Nrom(const double dtau, const Eigen::VectorXd& q, 
                        const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                        const Eigen::VectorXd& u);

  double residualSquaredNrom(const double dtau, const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                             const Eigen::VectorXd& u);

private:
  int time_step_, dimq_, dimv_;
  double barrier_, fraction_to_boundary_margin_;
  pdipmtest::JointVariablesUpperLimits position_upper_limits_, 
                                       velocity_upper_limits_, 
                                       torque_upper_limits_;
  pdipmtest::JointVariablesLowerLimits position_lower_limits_, 
                                       velocity_lower_limits_, 
                                       torque_lower_limits_;
};


inline JointSpaceConstraints::JointSpaceConstraints(const Robot& robot)
  : time_step_(10),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    barrier_(1.0e-08),
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


inline JointSpaceConstraints::JointSpaceConstraints()
  : time_step_(0),
    dimq_(0),
    dimv_(0),
    barrier_(0),
    fraction_to_boundary_margin_(0),
    position_upper_limits_(),
    velocity_upper_limits_(),
    torque_upper_limits_(),
    position_lower_limits_(),
    velocity_lower_limits_(),
    torque_lower_limits_() {
}


inline JointSpaceConstraints::~JointSpaceConstraints() {
}


inline void JointSpaceConstraints::setTimeStep(const int time_step) {
  assert(time_step >= 0);
  time_step_ = time_step;
}


inline bool JointSpaceConstraints::isFeasible(const Eigen::VectorXd& q, 
                                              const Eigen::VectorXd& v, 
                                              const Eigen::VectorXd& a, 
                                              const Eigen::VectorXd& u) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(u.size() == dimv_);
  if (!position_upper_limits_.isFeasible(q)) {
    return false;
  }
  if (!position_lower_limits_.isFeasible(q)) {
    return false;
  }
  if (!velocity_upper_limits_.isFeasible(v)) {
    return false;
  }
  if (!velocity_lower_limits_.isFeasible(v)) {
    return false;
  }
  if (!torque_upper_limits_.isFeasible(u)) {
    return false;
  }
  if (!torque_lower_limits_.isFeasible(u)) {
    return false;
  }
  return true;
}


inline void JointSpaceConstraints::setSlackAndDual(const double dtau, 
                                                   const Eigen::VectorXd& q, 
                                                   const Eigen::VectorXd& v, 
                                                   const Eigen::VectorXd& a, 
                                                   const Eigen::VectorXd& u) {
  assert(dtau > 0);
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(u.size() == dimv_);
  position_upper_limits_.setSlackAndDual(dtau, q);
  position_lower_limits_.setSlackAndDual(dtau, q);
  velocity_upper_limits_.setSlackAndDual(dtau, v);
  velocity_lower_limits_.setSlackAndDual(dtau, v);
  torque_upper_limits_.setSlackAndDual(dtau, u);
  torque_lower_limits_.setSlackAndDual(dtau, u);
}


inline void JointSpaceConstraints::augmentDualResidual(const double dtau, 
                                                       Eigen::VectorXd& Cu) {
  assert(dtau > 0);
  assert(Cu.size() == dimv_);
  torque_upper_limits_.augmentDualResidual(dtau, Cu);
  torque_lower_limits_.augmentDualResidual(dtau, Cu);
}


inline void JointSpaceConstraints::augmentDualResidual(const double dtau, 
                                                       Eigen::VectorXd& Cq,
                                                       Eigen::VectorXd& Cv, 
                                                       Eigen::VectorXd& Ca) {
  assert(dtau > 0);
  assert(Cq.size() == dimv_);
  assert(Cv.size() == dimv_);
  assert(Ca.size() == dimv_);
  position_upper_limits_.augmentDualResidual(dtau, Cq);
  position_lower_limits_.augmentDualResidual(dtau, Cq);
  velocity_upper_limits_.augmentDualResidual(dtau, Cv);
  velocity_lower_limits_.augmentDualResidual(dtau, Cv);
}


inline void JointSpaceConstraints::condenseSlackAndDual(const double dtau, 
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
  position_upper_limits_.condenseSlackAndDual(dtau, q, Cqq, Cq);
  position_lower_limits_.condenseSlackAndDual(dtau, q, Cqq, Cq);
  velocity_upper_limits_.condenseSlackAndDual(dtau, v, Cvv, Cv);
  velocity_lower_limits_.condenseSlackAndDual(dtau, v, Cvv, Cv);
}


inline void JointSpaceConstraints::condenseSlackAndDual(const double dtau, 
                                                        const Eigen::VectorXd& u, 
                                                        Eigen::MatrixXd& Cuu, 
                                                        Eigen::VectorXd& Cu) {
  assert(dtau > 0);
  assert(u.size() == dimv_);
  torque_upper_limits_.condenseSlackAndDual(dtau, u, Cuu, Cu);
  torque_lower_limits_.condenseSlackAndDual(dtau, u, Cuu, Cu);
}


inline void JointSpaceConstraints::computeSlackAndDualDirection(
    const double dtau, const Eigen::VectorXd& dq, const Eigen::VectorXd& dv, 
    const Eigen::VectorXd& da, const Eigen::VectorXd& du) {
  assert(dtau > 0);
  assert(dq.size() == dimv_);
  assert(dv.size() == dimv_);
  assert(da.size() == dimv_);
  assert(du.size() == dimv_);
  position_upper_limits_.computeSlackAndDualDirection(dtau, dq);
  position_lower_limits_.computeSlackAndDualDirection(dtau, dq);
  velocity_upper_limits_.computeSlackAndDualDirection(dtau, dv);
  velocity_lower_limits_.computeSlackAndDualDirection(dtau, dv);
  torque_upper_limits_.computeSlackAndDualDirection(dtau, du);
  torque_lower_limits_.computeSlackAndDualDirection(dtau, du);
}


inline double JointSpaceConstraints::maxSlackStepSize() {
  assert(fraction_to_boundary_margin_ > 0);
  assert(fraction_to_boundary_margin_ <= 1);
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


inline double JointSpaceConstraints::maxDualStepSize() {
  assert(fraction_to_boundary_margin_ > 0);
  assert(fraction_to_boundary_margin_ <= 1);
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


inline void JointSpaceConstraints::updateSlack(const double step_size) {
  assert(step_size > 0);
  position_upper_limits_.updateSlack(step_size);
  position_lower_limits_.updateSlack(step_size);
  velocity_upper_limits_.updateSlack(step_size);
  velocity_lower_limits_.updateSlack(step_size);
  torque_upper_limits_.updateSlack(step_size);
  torque_lower_limits_.updateSlack(step_size);
}


inline void JointSpaceConstraints::updateDual(const double step_size) {
  assert(step_size > 0);
  position_upper_limits_.updateDual(step_size);
  position_lower_limits_.updateDual(step_size);
  velocity_upper_limits_.updateDual(step_size);
  velocity_lower_limits_.updateDual(step_size);
  torque_upper_limits_.updateDual(step_size);
  torque_lower_limits_.updateDual(step_size);
}


inline double JointSpaceConstraints::costSlackBarrier() {
  double cost = 0;
  cost += position_upper_limits_.costSlackBarrier();
  cost += position_lower_limits_.costSlackBarrier();
  cost += velocity_upper_limits_.costSlackBarrier();
  cost += velocity_lower_limits_.costSlackBarrier();
  cost += torque_upper_limits_.costSlackBarrier();
  cost += torque_lower_limits_.costSlackBarrier();
  return cost;
}


inline double JointSpaceConstraints::costSlackBarrier(const double step_size) {
  double cost = 0;
  cost += position_upper_limits_.costSlackBarrier(step_size);
  cost += position_lower_limits_.costSlackBarrier(step_size);
  cost += velocity_upper_limits_.costSlackBarrier(step_size);
  cost += velocity_lower_limits_.costSlackBarrier(step_size);
  cost += torque_upper_limits_.costSlackBarrier(step_size);
  cost += torque_lower_limits_.costSlackBarrier(step_size);
  return cost;
}


inline double JointSpaceConstraints::residualL1Nrom(const double dtau, 
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
  norm += position_upper_limits_.residualL1Nrom(dtau, q);
  norm += position_lower_limits_.residualL1Nrom(dtau, q);
  norm += velocity_upper_limits_.residualL1Nrom(dtau, v);
  norm += velocity_lower_limits_.residualL1Nrom(dtau, v);
  norm += torque_upper_limits_.residualL1Nrom(dtau, u);
  norm += torque_lower_limits_.residualL1Nrom(dtau, u);
  return norm;
}


inline double JointSpaceConstraints::residualSquaredNrom(const double dtau, 
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
  norm += position_upper_limits_.residualSquaredNrom(dtau, q);
  norm += position_lower_limits_.residualSquaredNrom(dtau, q);
  norm += velocity_upper_limits_.residualSquaredNrom(dtau, v);
  norm += velocity_lower_limits_.residualSquaredNrom(dtau, v);
  norm += torque_upper_limits_.residualSquaredNrom(dtau, u);
  norm += torque_lower_limits_.residualSquaredNrom(dtau, u);
  return norm;
}


} // namespace pdipmtest
} // namespace idocp


#endif // IDOCP_CONSTRAINTS_JOINT_SPACE_CONSTRAINTS_HPP_