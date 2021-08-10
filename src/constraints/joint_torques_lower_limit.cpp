#include "idocp/constraints/joint_torques_lower_limit.hpp"


namespace idocp {

JointTorquesLowerLimit::JointTorquesLowerLimit(
    const Robot& robot, const double barrier, 
    const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
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


bool JointTorquesLowerLimit::useKinematics() const {
  return false;
}


KinematicsLevel JointTorquesLowerLimit::kinematicsLevel() const {
  return KinematicsLevel::AccelerationLevel;
}


bool JointTorquesLowerLimit::isFeasible(Robot& robot, 
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
                                      ConstraintComponentData& data, 
                                      const SplitSolution& s) const {
  data.slack = s.u - umin_;
}


void JointTorquesLowerLimit::computePrimalAndDualResidual(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s) const {
  data.residual = umin_ - s.u + data.slack;
  computeComplementarySlackness(data);
}


void JointTorquesLowerLimit::computePrimalResidualDerivatives(
    Robot& robot, ConstraintComponentData& data, const double dt, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  kkt_residual.lu.noalias() -= dt * data.dual;
}


void JointTorquesLowerLimit::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dt, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) const {
  kkt_matrix.Quu.diagonal().array()
      += dt * data.dual.array() / data.slack.array();
  kkt_residual.lu.array() 
      -= dt * (data.dual.array()*data.residual.array()-data.cmpl.array()) 
              / data.slack.array();
}


void JointTorquesLowerLimit::expandSlackAndDual(
    ConstraintComponentData& data, const SplitSolution& s, 
    const SplitDirection& d) const {
  data.dslack = d.du - data.residual;
  computeDualDirection(data);
}


int JointTorquesLowerLimit::dimc() const {
  return dimc_;
}

} // namespace idocp