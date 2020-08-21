#include "idocp/constraints/joint_acceleration_lower_limit.hpp"

#include <assert.h>


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


bool JointAccelerationLowerLimit::isFeasible(const Robot& robot, 
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
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  assert(dtau > 0);
  data.slack = dtau * (s.a.tail(dimc_)-amin_);
  setSlackAndDualPositive(data.slack, data.dual);
}


void JointAccelerationLowerLimit::augmentDualResidual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    KKTResidual& kkt_residual) const {
  kkt_residual.la().tail(dimc_).noalias() -= dtau * data.dual;
}


void JointAccelerationLowerLimit::condenseSlackAndDual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) const {
  kkt_matrix.Qaa().diagonal().tail(dimc_).array()
      += dtau * dtau * data.dual.array() / data.slack.array();
  data.residual = dtau * (amin_-s.a.tail(dimc_)) + data.slack;
  computeDuality(data.slack, data.dual, data.duality);
  kkt_residual.la().tail(dimc_).array() 
      -= dtau * (data.dual.array()*data.residual.array()-data.duality.array()) 
              / data.slack.array();
}


void JointAccelerationLowerLimit::computeSlackAndDualDirection(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitDirection& d) const {
  data.dslack = dtau * d.da().tail(dimc_) - data.residual;
  computeDualDirection(data.slack, data.dual, data.dslack, data.duality, 
                       data.ddual);
}


double JointAccelerationLowerLimit::residualL1Nrom(
    const Robot& robot, ConstraintComponentData& data, 
    const double dtau, const SplitSolution& s) const {
  data.residual = dtau * (amin_-s.a.tail(dimc_)) + data.slack;
  return data.residual.lpNorm<1>();
}


double JointAccelerationLowerLimit::squaredKKTErrorNorm(
    const Robot& robot, ConstraintComponentData& data, 
    const double dtau, const SplitSolution& s) const {
  data.residual = dtau * (amin_-s.a.tail(dimc_)) + data.slack;
  computeDuality(data.slack, data.dual, data.duality);
  double error = 0;
  error += data.residual.squaredNorm();
  error += data.duality.squaredNorm();
  return error;
}


int JointAccelerationLowerLimit::dimc() const {
  return dimc_;
}

} // namespace idocp