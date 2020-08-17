#include "idocp/constraints/joint_acceleration_upper_limit.hpp"

#include <assert.h>


namespace idocp {

JointAccelerationUpperLimit::JointAccelerationUpperLimit(
    const Robot& robot, const Eigen::VectorXd& amax, const double barrier, 
    const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(amax.size()),
    dim_passive_(robot.dim_passive()),
    amax_(amax) {
}


JointAccelerationUpperLimit::JointAccelerationUpperLimit()
  : ConstraintComponentBase(),
    dimc_(0),
    dim_passive_(0),
    amax_() {
}


JointAccelerationUpperLimit::~JointAccelerationUpperLimit() {
}


bool JointAccelerationUpperLimit::isFeasible(const Robot& robot, 
                                             ConstraintComponentData& data, 
                                             const SplitSolution& s) const {
  for (int i=0; i<dimc_; ++i) {
    if (s.a.tail(dimc_).coeff(i) > amax_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointAccelerationUpperLimit::setSlackAndDual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  assert(dtau > 0);
  data.slack = dtau * (amax_-s.a.tail(dimc_));
  setSlackAndDualPositive(data.slack, data.dual);
}


void JointAccelerationUpperLimit::augmentDualResidual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    KKTResidual& kkt_residual) const {
  kkt_residual.la().tail(dimc_).noalias() += dtau * data.dual;
}


void JointAccelerationUpperLimit::condenseSlackAndDual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) const {
  for (int i=0; i<dimc_; ++i) {
    kkt_matrix.Qaa().coeffRef(dim_passive_+i, dim_passive_+i) 
        += dtau * dtau * data.dual.coeff(i) / data.slack.coeff(i);
  }
  data.residual = dtau * (s.a.tail(dimc_)-amax_) + data.slack;
  computeDuality(data.slack, data.dual, data.duality);
  kkt_residual.la().tail(dimc_).array() 
      += dtau * (data.dual.array()*data.residual.array()-data.duality.array()) 
              / data.slack.array();
}


void JointAccelerationUpperLimit::computeSlackAndDualDirection(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitDirection& d) const {
  data.dslack = - dtau * d.da().tail(dimc_) - data.residual;
  computeDualDirection(data.slack, data.dual, data.dslack, data.duality, 
                       data.ddual);
}


double JointAccelerationUpperLimit::residualL1Nrom(
    const Robot& robot, ConstraintComponentData& data, 
    const double dtau, const SplitSolution& s) const {
  data.residual = dtau * (s.a.tail(dimc_)-amax_) + data.slack;
  return data.residual.lpNorm<1>();
}


double JointAccelerationUpperLimit::squaredKKTErrorNorm(
    const Robot& robot, ConstraintComponentData& data, 
    const double dtau, const SplitSolution& s) const {
  data.residual = dtau * (s.a.tail(dimc_)-amax_) + data.slack;
  computeDuality(data.slack, data.dual, data.duality);
  double error = 0;
  error += data.residual.squaredNorm();
  error += data.duality.squaredNorm();
  return error;
}


int JointAccelerationUpperLimit::dimc() const {
  return dimc_;
}

} // namespace idocp