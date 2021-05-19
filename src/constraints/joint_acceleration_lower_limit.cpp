#include "idocp/constraints/joint_acceleration_lower_limit.hpp"


namespace idocp {

JointAccelerationLowerLimit::JointAccelerationLowerLimit(
    const Robot& robot, const Eigen::VectorXd& amin, const double barrier, 
    const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(amin.size()),
    dim_passive_(robot.dim_passive()),
    amin_(amin) {
}


JointAccelerationLowerLimit::JointAccelerationLowerLimit()
  : ConstraintComponentBase(),
    dimc_(0),
    dim_passive_(0),
    amin_() {
}


JointAccelerationLowerLimit::~JointAccelerationLowerLimit() {
}


bool JointAccelerationLowerLimit::useKinematics() const {
  return false;
}


KinematicsLevel JointAccelerationLowerLimit::kinematicsLevel() const {
  return KinematicsLevel::AccelerationLevel;
}


bool JointAccelerationLowerLimit::isFeasible(Robot& robot, 
                                             ConstraintComponentData& data, 
                                             const SplitSolution& s) const {
  for (int i=0; i<dimc_; ++i) {
    if (s.a.tail(dimc_).coeff(i) < amin_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointAccelerationLowerLimit::setSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s) const {
  data.slack = s.a.tail(dimc_) - amin_;
  setSlackAndDualPositive(data);
}


void JointAccelerationLowerLimit::augmentDualResidual(
    Robot& robot, ConstraintComponentData& data, const double dt, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  kkt_residual.la.tail(dimc_).noalias() -= dt * data.dual;
}


void JointAccelerationLowerLimit::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dt, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) const {
  kkt_matrix.Qaa.diagonal().tail(dimc_).array()
      += dt * data.dual.array() / data.slack.array();
  computePrimalAndDualResidual(robot, data, s);
  kkt_residual.la.tail(dimc_).array() 
      -= dt * (data.dual.array()*data.residual.array()-data.duality.array()) 
              / data.slack.array();
}


void JointAccelerationLowerLimit::computeSlackAndDualDirection(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s, 
    const SplitDirection& d) const {
  data.dslack = d.da().tail(dimc_) - data.residual;
  computeDualDirection(data);
}


void JointAccelerationLowerLimit::computePrimalAndDualResidual(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s) const {
  data.residual = amin_ - s.a.tail(dimc_) + data.slack;
  computeDuality(data);
}


int JointAccelerationLowerLimit::dimc() const {
  return dimc_;
}

} // namespace idocp