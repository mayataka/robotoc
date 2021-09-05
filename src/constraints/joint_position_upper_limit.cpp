#include "idocp/constraints/joint_position_upper_limit.hpp"


namespace idocp {

JointPositionUpperLimit::JointPositionUpperLimit(
    const Robot& robot, const double barrier, 
    const double fraction_to_boundary_rule)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rule),
    dimc_(robot.lowerJointPositionLimit().size()),
    qmax_(robot.upperJointPositionLimit()) {
}


JointPositionUpperLimit::JointPositionUpperLimit()
  : ConstraintComponentBase(),
    dimc_(0),
    qmax_() {
}


JointPositionUpperLimit::~JointPositionUpperLimit() {
}


bool JointPositionUpperLimit::useKinematics() const {
  return false;
}

KinematicsLevel JointPositionUpperLimit::kinematicsLevel() const {
  return KinematicsLevel::PositionLevel;
}


bool JointPositionUpperLimit::isFeasible(Robot& robot, 
                                         ConstraintComponentData& data, 
                                         const SplitSolution& s) const {
  for (int i=0; i<dimc_; ++i) {
    if (s.q.tail(dimc_).coeff(i) > qmax_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointPositionUpperLimit::setSlack(Robot& robot, 
                                       ConstraintComponentData& data, 
                                       const SplitSolution& s) const {
  data.slack = qmax_ - s.q.tail(dimc_);
}


void JointPositionUpperLimit::evalConstraint(Robot& robot, 
                                             ConstraintComponentData& data, 
                                             const SplitSolution& s) const {
  data.residual = s.q.tail(dimc_) - qmax_ + data.slack;
  computeComplementarySlackness(data);
  data.log_barrier = logBarrier(data.slack);
}


void JointPositionUpperLimit::evalDerivatives(
    Robot& robot, ConstraintComponentData& data, const double dt, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  kkt_residual.lq().tail(dimc_).noalias() += dt * data.dual;
}


void JointPositionUpperLimit::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dt, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) const {
  kkt_matrix.Qqq().diagonal().tail(dimc_).array()
      += dt * data.dual.array() / data.slack.array();
  computeCondensingCoeffcient(data);
  kkt_residual.lq().tail(dimc_).noalias() += dt * data.cond;
}


void JointPositionUpperLimit::expandSlackAndDual(
    ConstraintComponentData& data, const SplitSolution& s, 
    const SplitDirection& d) const {
  data.dslack = - d.dq().tail(dimc_) - data.residual;
  computeDualDirection(data);
}


int JointPositionUpperLimit::dimc() const {
  return dimc_;
}

} // namespace idocp