#include "robotoc/constraints/joint_torques_upper_limit.hpp"


namespace robotoc {

JointTorquesUpperLimit::JointTorquesUpperLimit(const Robot& robot)
  : ConstraintComponentBase(),
    dimc_(robot.jointEffortLimit().size()),
    umax_(robot.jointEffortLimit()) {
}


JointTorquesUpperLimit::JointTorquesUpperLimit()
  : ConstraintComponentBase(),
    dimc_(0),
    umax_() {
}


JointTorquesUpperLimit::~JointTorquesUpperLimit() {
}


bool JointTorquesUpperLimit::useKinematics() const {
  return false;
}


KinematicsLevel JointTorquesUpperLimit::kinematicsLevel() const {
  return KinematicsLevel::AccelerationLevel;
}


bool JointTorquesUpperLimit::isFeasible(Robot& robot, 
                                        const ContactStatus& contact_status,
                                        ConstraintComponentData& data, 
                                        const SplitSolution& s) const {
  for (int i=0; i<dimc_; ++i) {
    if (s.u.coeff(i) > umax_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointTorquesUpperLimit::setSlack(Robot& robot, 
                                      const ContactStatus& contact_status,
                                      ConstraintComponentData& data, 
                                      const SplitSolution& s) const {
  data.slack = umax_ - s.u;
}


void JointTorquesUpperLimit::evalConstraint(Robot& robot, 
                                            const ContactStatus& contact_status,
                                            ConstraintComponentData& data, 
                                            const SplitSolution& s) const {
  data.residual = s.u - umax_ + data.slack;
  computeComplementarySlackness(data);
  data.log_barrier = logBarrier(data.slack);
}


void JointTorquesUpperLimit::evalDerivatives(
    Robot& robot, const ContactStatus& contact_status, 
    ConstraintComponentData& data, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) const {
  kkt_residual.lu.noalias() += data.dual;
}


void JointTorquesUpperLimit::condenseSlackAndDual( 
    const ContactStatus& contact_status, ConstraintComponentData& data, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) const {
  kkt_matrix.Quu.diagonal().array()
      += data.dual.array() / data.slack.array();
  computeCondensingCoeffcient(data);
  kkt_residual.lu.noalias() += data.cond;
}


void JointTorquesUpperLimit::expandSlackAndDual(
    const ContactStatus& contact_status, ConstraintComponentData& data, 
    const SplitDirection& d) const {
  data.dslack = - d.du - data.residual;
  computeDualDirection(data);
}


int JointTorquesUpperLimit::dimc() const {
  return dimc_;
}

} // namespace robotoc