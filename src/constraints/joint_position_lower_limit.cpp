#include "idocp/constraints/joint_position_lower_limit.hpp"

#include <assert.h>


namespace idocp {

JointPositionLowerLimit::JointPositionLowerLimit(
    const Robot& robot, const double barrier, 
    const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(robot.lowerJointPositionLimit().size()),
    dim_passive_(robot.dim_passive()),
    qmin_(robot.lowerJointPositionLimit()) {
}


JointPositionLowerLimit::JointPositionLowerLimit()
  : ConstraintComponentBase(),
    dimc_(0),
    dim_passive_(0),
    qmin_() {
}


JointPositionLowerLimit::~JointPositionLowerLimit() {
}


bool JointPositionLowerLimit::isFeasible(const Robot& robot, 
                                         ConstraintComponentData& data, 
                                         const SplitSolution& s) const {
  for (int i=0; i<dimc_; ++i) {
    if (s.q.tail(dimc_).coeff(i) < qmin_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointPositionLowerLimit::setSlackAndDual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  assert(dtau > 0);
  data.slack = dtau * (s.q.tail(dimc_)-qmin_);
  setSlackAndDualPositive(data.slack, data.dual);
}


void JointPositionLowerLimit::augmentDualResidual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    KKTResidual& kkt_residual) const {
  kkt_residual.lq().tail(dimc_).noalias() -= dtau * data.dual;
}


void JointPositionLowerLimit::condenseSlackAndDual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) const {
  for (int i=0; i<dimc_; ++i) {
    kkt_matrix.Qqq().coeffRef(dim_passive_+i, dim_passive_+i) 
        += dtau * dtau * data.dual.coeff(i) / data.slack.coeff(i);
  }
  data.residual = dtau * (qmin_-s.q.tail(dimc_)) + data.slack;
  computeDuality(data.slack, data.dual, data.duality);
  kkt_residual.lq().tail(dimc_).array() 
      -= dtau * (data.dual.array()*data.residual.array()-data.duality.array()) 
              / data.slack.array();
}


void JointPositionLowerLimit::computeSlackAndDualDirection(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitDirection& d) const {
  data.dslack = dtau * d.dq().tail(dimc_) - data.residual;
  computeDualDirection(data.slack, data.dual, data.dslack, data.duality, 
                       data.ddual);
}


double JointPositionLowerLimit::residualL1Nrom(
    const Robot& robot, ConstraintComponentData& data, 
    const double dtau, const SplitSolution& s) const {
  data.residual = dtau * (qmin_-s.q.tail(dimc_)) + data.slack;
  return data.residual.lpNorm<1>();
}


double JointPositionLowerLimit::squaredKKTErrorNorm(
    const Robot& robot, ConstraintComponentData& data, 
    const double dtau, const SplitSolution& s) const {
  data.residual = dtau * (qmin_-s.q.tail(dimc_)) + data.slack;
  computeDuality(data.slack, data.dual, data.duality);
  double error = 0;
  error += data.residual.squaredNorm();
  error += data.duality.squaredNorm();
  return error;
}


int JointPositionLowerLimit::dimc() const {
  return dimc_;
}

} // namespace idocp