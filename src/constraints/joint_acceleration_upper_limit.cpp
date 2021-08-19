#include "idocp/constraints/joint_acceleration_upper_limit.hpp"


namespace idocp {

JointAccelerationUpperLimit::JointAccelerationUpperLimit(
    const Robot& robot, const Eigen::VectorXd& amax, const double barrier, 
    const double fraction_to_boundary_rule)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rule),
    dimc_(amax.size()),
    amax_(amax) {
}


JointAccelerationUpperLimit::JointAccelerationUpperLimit()
  : ConstraintComponentBase(),
    dimc_(0),
    amax_() {
}


JointAccelerationUpperLimit::~JointAccelerationUpperLimit() {
}


bool JointAccelerationUpperLimit::useKinematics() const {
  return false;
}


KinematicsLevel JointAccelerationUpperLimit::kinematicsLevel() const {
  return KinematicsLevel::AccelerationLevel;
}


bool JointAccelerationUpperLimit::isFeasible(Robot& robot, 
                                             ConstraintComponentData& data, 
                                             const SplitSolution& s) const {
  for (int i=0; i<dimc_; ++i) {
    if (s.a.tail(dimc_).coeff(i) > amax_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointAccelerationUpperLimit::setSlack(Robot& robot, 
                                           ConstraintComponentData& data, 
                                           const SplitSolution& s) const {
  data.slack = amax_ - s.a.tail(dimc_);
}


void JointAccelerationUpperLimit::computePrimalAndDualResidual(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s) const {
  data.residual = s.a.tail(dimc_) - amax_ + data.slack;
  computeComplementarySlackness(data);
  data.log_barrier = logBarrier(data.slack);
}


void JointAccelerationUpperLimit::computePrimalResidualDerivatives(
    Robot& robot, ConstraintComponentData& data, const double dt, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  kkt_residual.la.tail(dimc_).noalias() += dt * data.dual;
}


void JointAccelerationUpperLimit::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dt, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) const {
  kkt_matrix.Qaa.diagonal().tail(dimc_).array()
      += dt * data.dual.array() / data.slack.array();
  computeCondensingCoeffcient(data);
  kkt_residual.la.tail(dimc_).noalias() += dt * data.cond;
}


void JointAccelerationUpperLimit::expandSlackAndDual(
    ConstraintComponentData& data, const SplitSolution& s, 
    const SplitDirection& d) const {
  data.dslack = - d.da().tail(dimc_) - data.residual;
  computeDualDirection(data);
}


int JointAccelerationUpperLimit::dimc() const {
  return dimc_;
}

} // namespace idocp