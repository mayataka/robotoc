#include "idocp/constraints/joint_torques_upper_limit.hpp"

#include <assert.h>


namespace idocp {

JointTorquesUpperLimit::JointTorquesUpperLimit(
    const Robot& robot, const double barrier, 
    const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(robot.jointEffortLimit().size()),
    dim_passive_(robot.dim_passive()),
    umax_(robot.jointEffortLimit()) {
}


JointTorquesUpperLimit::JointTorquesUpperLimit()
  : ConstraintComponentBase(),
    dimc_(0),
    dim_passive_(0),
    umax_() {
}


JointTorquesUpperLimit::~JointTorquesUpperLimit() {
}


bool JointTorquesUpperLimit::isFeasible(const Robot& robot, 
                                        ConstraintComponentData& data, 
                                        const SplitSolution& s) const {
  for (int i=0; i<dimc_; ++i) {
    if (s.u.tail(dimc_).coeff(i) > umax_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointTorquesUpperLimit::setSlackAndDual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  assert(dtau > 0);
  data.slack = dtau * (umax_-s.u.tail(dimc_));
  setSlackAndDualPositive(data.slack, data.dual);
}


void JointTorquesUpperLimit::augmentDualResidual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    KKTResidual& kkt_residual) const {
  kkt_residual.lu.tail(dimc_).noalias() += dtau * data.dual;
}


void JointTorquesUpperLimit::condenseSlackAndDual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) const {
  for (int i=0; i<dimc_; ++i) {
    kkt_matrix.Quu.coeffRef(dim_passive_+i, dim_passive_+i) 
        += dtau * dtau * data.dual.coeff(i) / data.slack.coeff(i);
  }
  data.residual = dtau * (s.u.tail(dimc_)-umax_) + data.slack;
  computeDuality(data.slack, data.dual, data.duality);
  kkt_residual.lu.tail(dimc_).array() 
      += dtau * (data.dual.array()*data.residual.array()-data.duality.array()) 
              / data.slack.array();
}


void JointTorquesUpperLimit::computeSlackAndDualDirection(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitDirection& d) const {
  data.dslack = - dtau * d.du.tail(dimc_) - data.residual;
  computeDualDirection(data.slack, data.dual, data.dslack, data.duality, 
                       data.ddual);
}


double JointTorquesUpperLimit::residualL1Nrom(
    const Robot& robot, ConstraintComponentData& data, 
    const double dtau, const SplitSolution& s) const {
  data.residual = dtau * (s.u.tail(dimc_)-umax_) + data.slack;
  return data.residual.lpNorm<1>();
}


double JointTorquesUpperLimit::squaredKKTErrorNorm(
    const Robot& robot, ConstraintComponentData& data, 
    const double dtau, const SplitSolution& s) const {
  data.residual = dtau * (s.u.tail(dimc_)-umax_) + data.slack;
  computeDuality(data.slack, data.dual, data.duality);
  double error = 0;
  error += data.residual.squaredNorm();
  error += data.duality.squaredNorm();
  return error;
}


int JointTorquesUpperLimit::dimc() const {
  return dimc_;
}

} // namespace idocp