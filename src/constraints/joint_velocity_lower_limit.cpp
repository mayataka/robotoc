#include "idocp/constraints/joint_velocity_lower_limit.hpp"


namespace idocp {

JointVelocityLowerLimit::JointVelocityLowerLimit(
    const Robot& robot, const double barrier, 
    const double fraction_to_boundary_rule)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rule),
    dimc_(robot.jointVelocityLimit().size()),
    vmin_(-robot.jointVelocityLimit()) {
}


JointVelocityLowerLimit::JointVelocityLowerLimit()
  : ConstraintComponentBase(),
    dimc_(0),
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


void JointVelocityLowerLimit::setSlack(Robot& robot, 
                                       ConstraintComponentData& data, 
                                       const SplitSolution& s) const {
  data.slack = s.v.tail(dimc_) - vmin_;
}


void JointVelocityLowerLimit::computePrimalAndDualResidual(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s) const {
  data.residual = vmin_ - s.v.tail(dimc_) + data.slack;
  computeComplementarySlackness(data);
  data.log_barrier = logBarrier(data.slack);
}


void JointVelocityLowerLimit::computePrimalResidualDerivatives(
    Robot& robot, ConstraintComponentData& data, const double dt, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  kkt_residual.lv().tail(dimc_).noalias() -= dt * data.dual;
}


void JointVelocityLowerLimit::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dt, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) const {
  kkt_matrix.Qvv().diagonal().tail(dimc_).array()
      += dt * data.dual.array() / data.slack.array();
  kkt_residual.lv().tail(dimc_).array() 
      -= dt * (data.dual.array()*data.residual.array()-data.cmpl.array()) 
              / data.slack.array();
}


void JointVelocityLowerLimit::expandSlackAndDual(
    ConstraintComponentData& data, const SplitSolution& s, 
    const SplitDirection& d) const {
  data.dslack = d.dv().tail(dimc_) - data.residual;
  computeDualDirection(data);
}


int JointVelocityLowerLimit::dimc() const {
  return dimc_;
}

} // namespace idocp