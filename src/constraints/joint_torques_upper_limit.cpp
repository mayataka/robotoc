#include "robotoc/constraints/joint_torques_upper_limit.hpp"


namespace robotoc {

JointTorquesUpperLimit::JointTorquesUpperLimit(
    const Robot& robot, const double barrier, 
    const double fraction_to_boundary_rule)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rule),
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
                                      ConstraintComponentData& data, 
                                      const SplitSolution& s) const {
  data.slack = umax_ - s.u;
}


void JointTorquesUpperLimit::evalConstraint(Robot& robot, 
                                                ConstraintComponentData& data, 
                                                const SplitSolution& s) const {
  data.residual = s.u - umax_ + data.slack;
  computeComplementarySlackness(data);
  data.log_barrier = logBarrier(data.slack);
}


void JointTorquesUpperLimit::evalDerivatives(
    Robot& robot, ConstraintComponentData& data, const double dt, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  kkt_residual.lu.noalias() += dt * data.dual;
}


void JointTorquesUpperLimit::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dt, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) const {
  kkt_matrix.Quu.diagonal().array()
      += dt * data.dual.array() / data.slack.array();
  computeCondensingCoeffcient(data);
  kkt_residual.lu.noalias() += dt * data.cond;
}


void JointTorquesUpperLimit::expandSlackAndDual(
    ConstraintComponentData& data, const SplitSolution& s, 
    const SplitDirection& d) const {
  data.dslack = - d.du - data.residual;
  computeDualDirection(data);
}


int JointTorquesUpperLimit::dimc() const {
  return dimc_;
}

} // namespace robotoc