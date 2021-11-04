#include "robotoc/constraints/joint_position_lower_limit.hpp"


namespace robotoc {

JointPositionLowerLimit::JointPositionLowerLimit(
    const Robot& robot, const double barrier, 
    const double fraction_to_boundary_rule)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rule),
    dimc_(robot.lowerJointPositionLimit().size()),
    qmin_(robot.lowerJointPositionLimit()) {
}


JointPositionLowerLimit::JointPositionLowerLimit()
  : ConstraintComponentBase(),
    dimc_(0),
    qmin_() {
}


JointPositionLowerLimit::~JointPositionLowerLimit() {
}


bool JointPositionLowerLimit::useKinematics() const {
  return false;
}


KinematicsLevel JointPositionLowerLimit::kinematicsLevel() const {
  return KinematicsLevel::PositionLevel;
}


bool JointPositionLowerLimit::isFeasible(Robot& robot, 
                                         const ContactStatus& contact_status,
                                         ConstraintComponentData& data, 
                                         const SplitSolution& s) const {
  for (int i=0; i<dimc_; ++i) {
    if (s.q.tail(dimc_).coeff(i) < qmin_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointPositionLowerLimit::setSlack(Robot& robot, 
                                       const ContactStatus& contact_status,
                                       ConstraintComponentData& data, 
                                       const SplitSolution& s) const {
  data.slack = s.q.tail(dimc_) - qmin_;
}


void JointPositionLowerLimit::evalConstraint(Robot& robot, 
                                             const ContactStatus& contact_status,
                                             ConstraintComponentData& data, 
                                             const SplitSolution& s) const {
  data.residual = qmin_ - s.q.tail(dimc_) + data.slack;
  computeComplementarySlackness(data);
  data.log_barrier = logBarrier(data.slack);
}


void JointPositionLowerLimit::evalDerivatives(
    Robot& robot, const ContactStatus& contact_status,
    ConstraintComponentData& data, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) const {
  kkt_residual.lq().tail(dimc_).noalias() -= data.dual;
}


void JointPositionLowerLimit::condenseSlackAndDual(
    const ContactStatus& contact_status, ConstraintComponentData& data, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) const {
  kkt_matrix.Qqq().diagonal().tail(dimc_).array()
      += data.dual.array() / data.slack.array();
  computeCondensingCoeffcient(data);
  kkt_residual.lq().tail(dimc_).noalias() -= data.cond;
}


void JointPositionLowerLimit::expandSlackAndDual(
    const ContactStatus& contact_status, ConstraintComponentData& data, 
    const SplitDirection& d) const {
  data.dslack = d.dq().tail(dimc_) - data.residual;
  computeDualDirection(data);
}


int JointPositionLowerLimit::dimc() const {
  return dimc_;
}

} // namespace robotoc