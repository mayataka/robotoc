#include "idocp/constraints/joint_velocity_lower_limit.hpp"

#include <assert.h>


namespace idocp {

JointVelocityLowerLimit::JointVelocityLowerLimit(
    const Robot& robot, const double barrier, 
    const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(robot.jointVelocityLimit().size()),
    dim_passive_(robot.dim_passive()),
    vmin_(-robot.jointVelocityLimit()) {
}


JointVelocityLowerLimit::JointVelocityLowerLimit()
  : ConstraintComponentBase(),
    dimc_(0),
    dim_passive_(0),
    vmin_() {
}


JointVelocityLowerLimit::~JointVelocityLowerLimit() {
}


bool JointVelocityLowerLimit::useKinematics() const {
  return false;
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


void JointVelocityLowerLimit::setSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  assert(dtau > 0);
  data.slack = dtau * (s.v.tail(dimc_)-vmin_);
  setSlackAndDualPositive(data);
}


void JointVelocityLowerLimit::augmentDualResidual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, KKTResidual& kkt_residual) const {
  kkt_residual.lv().tail(dimc_).noalias() -= dtau * data.dual;
}


void JointVelocityLowerLimit::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) const {
  kkt_matrix.Qvv().diagonal().tail(dimc_).array()
      += dtau * dtau * data.dual.array() / data.slack.array();
  data.residual = dtau * (vmin_-s.v.tail(dimc_)) + data.slack;
  computeDuality(data);
  kkt_residual.lv().tail(dimc_).array() 
      -= dtau * (data.dual.array()*data.residual.array()-data.duality.array()) 
              / data.slack.array();
}


void JointVelocityLowerLimit::computeSlackAndDualDirection(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, const SplitDirection& d) const {
  data.dslack = dtau * d.dv().tail(dimc_) - data.residual;
  computeDualDirection(data);
}


double JointVelocityLowerLimit::residualL1Nrom(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  data.residual = dtau * (vmin_-s.v.tail(dimc_)) + data.slack;
  return data.residual.lpNorm<1>();
}


double JointVelocityLowerLimit::squaredKKTErrorNorm(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  data.residual = dtau * (vmin_-s.v.tail(dimc_)) + data.slack;
  computeDuality(data);
  double error = 0;
  error += data.residual.squaredNorm();
  error += data.duality.squaredNorm();
  return error;
}


int JointVelocityLowerLimit::dimc() const {
  return dimc_;
}

} // namespace idocp