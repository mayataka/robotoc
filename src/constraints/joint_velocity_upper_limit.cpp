#include "idocp/constraints/joint_velocity_upper_limit.hpp"

#include <assert.h>


namespace idocp {

JointVelocityUpperLimit::JointVelocityUpperLimit(
    const Robot& robot, const double barrier, 
    const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(robot.jointVelocityLimit().size()),
    dim_passive_(robot.dim_passive()),
    vmax_(robot.jointVelocityLimit()) {
}


JointVelocityUpperLimit::JointVelocityUpperLimit()
  : ConstraintComponentBase(),
    dimc_(0),
    dim_passive_(0),
    vmax_() {
}


JointVelocityUpperLimit::~JointVelocityUpperLimit() {
}


bool JointVelocityUpperLimit::useKinematics() const {
  return false;
}


bool JointVelocityUpperLimit::isFeasible(Robot& robot, 
                                         ConstraintComponentData& data, 
                                         const SplitSolution& s) const {
  for (int i=0; i<dimc_; ++i) {
    if (s.v.tail(dimc_).coeff(i) > vmax_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointVelocityUpperLimit::setSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  assert(dtau > 0);
  data.slack = dtau * (vmax_-s.v.tail(dimc_));
  setSlackAndDualPositive(data.slack, data.dual);
}


void JointVelocityUpperLimit::augmentDualResidual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    KKTResidual& kkt_residual) const {
  kkt_residual.lv().tail(dimc_).noalias() += dtau * data.dual;
}


void JointVelocityUpperLimit::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) const {
  kkt_matrix.Qvv().diagonal().tail(dimc_).array()
      += dtau * dtau * data.dual.array() / data.slack.array();
  data.residual = dtau * (s.v.tail(dimc_)-vmax_) + data.slack;
  computeDuality(data.slack, data.dual, data.duality);
  kkt_residual.lv().tail(dimc_).array() 
      += dtau * (data.dual.array()*data.residual.array()-data.duality.array()) 
              / data.slack.array();
}


void JointVelocityUpperLimit::computeSlackAndDualDirection(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitDirection& d) const {
  data.dslack = - dtau * d.dv().tail(dimc_) - data.residual;
  computeDualDirection(data.slack, data.dual, data.dslack, data.duality, 
                       data.ddual);
}


double JointVelocityUpperLimit::residualL1Nrom(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  data.residual = dtau * (s.v.tail(dimc_)-vmax_) + data.slack;
  return data.residual.lpNorm<1>();
}


double JointVelocityUpperLimit::squaredKKTErrorNorm(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  data.residual = dtau * (s.v.tail(dimc_)-vmax_) + data.slack;
  computeDuality(data.slack, data.dual, data.duality);
  double error = 0;
  error += data.residual.squaredNorm();
  error += data.duality.squaredNorm();
  return error;
}


int JointVelocityUpperLimit::dimc() const {
  return dimc_;
}

} // namespace idocp