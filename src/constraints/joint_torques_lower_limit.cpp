#include "idocp/constraints/joint_torques_lower_limit.hpp"

#include <assert.h>


namespace idocp {

JointTorquesLowerLimit::JointTorquesLowerLimit(
    const Robot& robot, const double barrier, 
    const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(robot.jointEffortLimit().size()),
    dim_passive_(robot.dim_passive()),
    umin_(-robot.jointEffortLimit()) {
}


JointTorquesLowerLimit::JointTorquesLowerLimit()
  : ConstraintComponentBase(),
    dimc_(0),
    dim_passive_(0),
    umin_() {
}


JointTorquesLowerLimit::~JointTorquesLowerLimit() {
}


bool JointTorquesLowerLimit::isFeasible(const Robot& robot, 
                                        ConstraintComponentData& data, 
                                        const SplitSolution& s) const {
  for (int i=0; i<dimc_; ++i) {
    if (s.u.tail(dimc_).coeff(i) < umin_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointTorquesLowerLimit::setSlackAndDual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  assert(dtau > 0);
  data.slack = dtau * (s.u.tail(dimc_)-umin_);
  setSlackAndDualPositive(data.slack, data.dual);
}


void JointTorquesLowerLimit::augmentDualResidual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    Eigen::VectorXd& lu) const {
  lu.tail(dimc_).noalias() -= dtau * data.dual;
}


void JointTorquesLowerLimit::condenseSlackAndDual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const Eigen::VectorXd& u, Eigen::MatrixXd& Quu, Eigen::VectorXd& lu) const {
  Quu.diagonal().tail(dimc_).array()
      += dtau * dtau * data.dual.array() / data.slack.array();
  data.residual = dtau * (umin_-u.tail(dimc_)) + data.slack;
  computeDuality(data.slack, data.dual, data.duality);
  lu.tail(dimc_).array() 
      -= dtau * (data.dual.array()*data.residual.array()-data.duality.array()) 
              / data.slack.array();
}


void JointTorquesLowerLimit::computeSlackAndDualDirection(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitDirection& d) const {
  data.dslack = dtau * d.du.tail(dimc_) - data.residual;
  computeDualDirection(data.slack, data.dual, data.dslack, data.duality, 
                       data.ddual);
}


double JointTorquesLowerLimit::residualL1Nrom(
    const Robot& robot, ConstraintComponentData& data, 
    const double dtau, const SplitSolution& s) const {
  data.residual = dtau * (umin_-s.u.tail(dimc_)) + data.slack;
  return data.residual.lpNorm<1>();
}


double JointTorquesLowerLimit::squaredKKTErrorNorm(
    const Robot& robot, ConstraintComponentData& data, 
    const double dtau, const SplitSolution& s) const {
  data.residual = dtau * (umin_-s.u.tail(dimc_)) + data.slack;
  computeDuality(data.slack, data.dual, data.duality);
  double error = 0;
  error += data.residual.squaredNorm();
  error += data.duality.squaredNorm();
  return error;
}


int JointTorquesLowerLimit::dimc() const {
  return dimc_;
}

} // namespace idocp