#include "idocp/constraints/joint_velocity_lower_limit.hpp"

#include <cassert>


namespace idocp {

JointVelocityLowerLimit::JointVelocityLowerLimit(
    const Robot& robot, const double barrier, 
    const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(robot.jointVelocityLimit().size()),
    dim_passive_(robot.dim_passive()),
    vmin_(-robot.jointVelocityLimit()) {
}


JointVelocityLowerLimit::JointVelocityLowerLimit()
  : ConstraintComponentBase(),
    dimc_(0),
    dim_passive_(0),
    vmin_() {
}


JointVelocityLowerLimit::~JointVelocityLowerLimit() {
}


bool JointVelocityLowerLimit::useKinematics() const {
  return false;
}


KinematicsLevel JointVelocityLowerLimit::kinematicsLevel() const {
  return KinematicsLevel::VelocityLevel;
}


bool JointVelocityLowerLimit::isFeasible(Robot& robot, 
                                         ConstraintComponentData& data, 
                                         const SplitSolution& s) const {
  for (int i=0; i<dimc_; ++i) {
    if (s.v.tail(dimc_).coeff(i) < vmin_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointVelocityLowerLimit::setSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s) const {
  data.slack = s.v.tail(dimc_) - vmin_;
  setSlackAndDualPositive(data);
}


void JointVelocityLowerLimit::augmentDualResidual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  assert(dtau > 0);
  kkt_residual.lv().tail(dimc_).noalias() -= dtau * data.dual;
}


void JointVelocityLowerLimit::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) const {
  assert(dtau > 0);
  kkt_matrix.Qvv().diagonal().tail(dimc_).array()
      += dtau * data.dual.array() / data.slack.array();
  computePrimalAndDualResidual(robot, data, s);
  kkt_residual.lv().tail(dimc_).array() 
      -= dtau * (data.dual.array()*data.residual.array()-data.duality.array()) 
              / data.slack.array();
}


void JointVelocityLowerLimit::computeSlackAndDualDirection(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s, 
    const SplitDirection& d) const {
  data.dslack = d.dv().tail(dimc_) - data.residual;
  computeDualDirection(data);
}


void JointVelocityLowerLimit::computePrimalAndDualResidual(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s) const {
  data.residual = vmin_ - s.v.tail(dimc_) + data.slack;
  computeDuality(data);
}


int JointVelocityLowerLimit::dimc() const {
  return dimc_;
}

} // namespace idocp