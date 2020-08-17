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


bool JointVelocityLowerLimit::isFeasible(const Robot& robot, 
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
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  assert(dtau > 0);
  data.slack = dtau * (s.v.tail(dimc_)-vmin_);
  setSlackAndDualPositive(data.slack, data.dual);
}


void JointVelocityLowerLimit::augmentDualResidual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    KKTResidual& kkt_residual) const {
  kkt_residual.lv().tail(dimc_).noalias() -= dtau * data.dual;
}


void JointVelocityLowerLimit::condenseSlackAndDual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) const {
  for (int i=0; i<dimc_; ++i) {
    kkt_matrix.Qvv().coeffRef(dim_passive_+i, dim_passive_+i) 
        += dtau * dtau * data.dual.coeff(i) / data.slack.coeff(i);
  }
  data.residual = dtau * (vmin_-s.v.tail(dimc_)) + data.slack;
  computeDualityResidual(data.slack, data.dual, data.duality);
  kkt_residual.lv().tail(dimc_).array() 
      -= dtau * (data.dual.array()*data.residual.array()-data.duality.array()) 
              / data.slack.array();
}


void JointVelocityLowerLimit::computeSlackAndDualDirection(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitDirection& d) const {
  data.dslack = dtau * d.dv().tail(dimc_) - data.residual;
  computeDualDirection(data.slack, data.dslack, data.dual, data.duality, 
                       data.ddual);
}


double JointVelocityLowerLimit::residualL1Nrom(
    const Robot& robot, ConstraintComponentData& data, 
    const double dtau, const SplitSolution& s) const {
  data.residual = dtau * (vmin_-s.v.tail(dimc_)) + data.slack;
  return data.residual.lpNorm<1>();
}


double JointVelocityLowerLimit::squaredKKTErrorNorm(
    const Robot& robot, ConstraintComponentData& data, 
    const double dtau, const SplitSolution& s) const {
  data.residual = dtau * (vmin_-s.v.tail(dimc_)) + data.slack;
  computeDualityResidual(data.slack, data.dual, data.duality);
  double error = 0;
  error += data.residual.squaredNorm();
  error += data.duality.squaredNorm();
  return error;
}


int JointVelocityLowerLimit::dimc() const {
  return dimc_;
}

} // namespace idocp