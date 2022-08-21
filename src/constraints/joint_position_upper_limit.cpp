#include "robotoc/constraints/joint_position_upper_limit.hpp"


namespace robotoc {

JointPositionUpperLimit::JointPositionUpperLimit(const Robot& robot)
  : ConstraintComponentBase(),
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


KinematicsLevel JointPositionUpperLimit::kinematicsLevel() const {
  return KinematicsLevel::PositionLevel;
}


bool JointPositionUpperLimit::isFeasible(Robot& robot, 
                                         const ContactStatus& contact_status,
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
                                       const ContactStatus& contact_status,
                                       ConstraintComponentData& data, 
                                       const SplitSolution& s) const {
  data.slack = qmax_ - s.q.tail(dimc_);
}


void JointPositionUpperLimit::evalConstraint(Robot& robot, 
                                             const ContactStatus& contact_status,
                                             ConstraintComponentData& data, 
                                             const SplitSolution& s) const {
  data.residual = s.q.tail(dimc_) - qmax_ + data.slack;
  computeComplementarySlackness(data);
  data.log_barrier = logBarrier(data.slack);
}


void JointPositionUpperLimit::evalDerivatives(
    Robot& robot, const ContactStatus& contact_status, 
    ConstraintComponentData& data, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) const {
  kkt_residual.lq().tail(dimc_).noalias() += data.dual;
}


void JointPositionUpperLimit::condenseSlackAndDual(
    const ContactStatus& contact_status, ConstraintComponentData& data, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) const {
  kkt_matrix.Qqq().diagonal().tail(dimc_).array()
      += data.dual.array() / data.slack.array();
  computeCondensingCoeffcient(data);
  kkt_residual.lq().tail(dimc_).noalias() += data.cond;
}


void JointPositionUpperLimit::expandSlackAndDual(
    const ContactStatus& contact_status, ConstraintComponentData& data, 
    const SplitDirection& d) const {
  data.dslack = - d.dq().tail(dimc_) - data.residual;
  computeDualDirection(data);
}


int JointPositionUpperLimit::dimc() const {
  return dimc_;
}

} // namespace robotoc