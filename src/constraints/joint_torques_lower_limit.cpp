#include "robotoc/constraints/joint_torques_lower_limit.hpp"


namespace robotoc {

JointTorquesLowerLimit::JointTorquesLowerLimit(const Robot& robot)
  : ConstraintComponentBase(),
    dimc_(robot.jointEffortLimit().size()),
    umin_(-robot.jointEffortLimit()) {
}


JointTorquesLowerLimit::JointTorquesLowerLimit()
  : ConstraintComponentBase(),
    dimc_(0),
    umin_() {
}


JointTorquesLowerLimit::~JointTorquesLowerLimit() {
}


KinematicsLevel JointTorquesLowerLimit::kinematicsLevel() const {
  return KinematicsLevel::AccelerationLevel;
}


bool JointTorquesLowerLimit::isFeasible(Robot& robot, 
                                        const ContactStatus& contact_status, 
                                        ConstraintComponentData& data, 
                                        const SplitSolution& s) const {
  for (int i=0; i<dimc_; ++i) {
    if (s.u.coeff(i) < umin_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointTorquesLowerLimit::setSlack(Robot& robot, 
                                      const ContactStatus& contact_status, 
                                      ConstraintComponentData& data, 
                                      const SplitSolution& s) const {
  data.slack = s.u - umin_;
}


void JointTorquesLowerLimit::evalConstraint(Robot& robot, 
                                            const ContactStatus& contact_status, 
                                            ConstraintComponentData& data, 
                                            const SplitSolution& s) const {
  data.residual = umin_ - s.u + data.slack;
  computeComplementarySlackness(data);
  data.log_barrier = logBarrier(data.slack);
}


void JointTorquesLowerLimit::evalDerivatives(
    Robot& robot, const ContactStatus& contact_status, 
    ConstraintComponentData& data, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) const {
  kkt_residual.lu.noalias() -= data.dual;
}


void JointTorquesLowerLimit::condenseSlackAndDual(
    const ContactStatus& contact_status, ConstraintComponentData& data, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) const {
  kkt_matrix.Quu.diagonal().array()
      += data.dual.array() / data.slack.array();
  computeCondensingCoeffcient(data);
  kkt_residual.lu.noalias() -= data.cond;
}


void JointTorquesLowerLimit::expandSlackAndDual(
    const ContactStatus& contact_status, ConstraintComponentData& data, 
    const SplitDirection& d) const {
  data.dslack = d.du - data.residual;
  computeDualDirection(data);
}


int JointTorquesLowerLimit::dimc() const {
  return dimc_;
}

} // namespace robotoc