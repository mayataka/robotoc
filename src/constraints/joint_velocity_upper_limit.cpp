#include "idocp/constraints/joint_velocity_upper_limit.hpp"

#include <cassert>


namespace idocp {

JointVelocityUpperLimit::JointVelocityUpperLimit(
    const Robot& robot, const double barrier, 
    const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(robot.jointVelocityLimit().size()),
    dim_passive_(robot.dim_passive()),
    vmax_(robot.jointVelocityLimit()) {
}


JointVelocityUpperLimit::JointVelocityUpperLimit()
  : ConstraintComponentBase(),
    dimc_(0),
    dim_passive_(0),
    vmax_() {
}


JointVelocityUpperLimit::~JointVelocityUpperLimit() {
}


bool JointVelocityUpperLimit::useKinematics() const {
  return false;
}


KinematicsLevel JointVelocityUpperLimit::kinematicsLevel() const {
  return KinematicsLevel::VelocityLevel;
}


bool JointVelocityUpperLimit::isFeasible(Robot& robot, 
                                         ConstraintComponentData& data, 
                                         const SplitSolution& s) const {
  for (int i=0; i<dimc_; ++i) {
    if (s.v.tail(dimc_).coeff(i) > vmax_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointVelocityUpperLimit::setSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s) const {
  data.slack = vmax_ - s.v.tail(dimc_);
  setSlackAndDualPositive(data);
}


void JointVelocityUpperLimit::augmentDualResidual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  assert(dtau > 0);
  kkt_residual.lv().tail(dimc_).noalias() += dtau * data.dual;
}


void JointVelocityUpperLimit::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) const {
  assert(dtau > 0);
  kkt_matrix.Qvv().diagonal().tail(dimc_).array()
      += dtau * data.dual.array() / data.slack.array();
  computePrimalAndDualResidual(robot, data, s);
  kkt_residual.lv().tail(dimc_).array() 
      += dtau * (data.dual.array()*data.residual.array()-data.duality.array()) 
              / data.slack.array();
}


void JointVelocityUpperLimit::computeSlackAndDualDirection(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s, 
    const SplitDirection& d) const {
  data.dslack = - d.dv().tail(dimc_) - data.residual;
  computeDualDirection(data);
}


void JointVelocityUpperLimit::computePrimalAndDualResidual(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s) const {
  data.residual = s.v.tail(dimc_) - vmax_ + data.slack;
  computeDuality(data);
}


int JointVelocityUpperLimit::dimc() const {
  return dimc_;
}

} // namespace idocp