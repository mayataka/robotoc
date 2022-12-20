#include "robotoc/constraints/joint_acceleration_upper_limit.hpp"


namespace robotoc {

JointAccelerationUpperLimit::JointAccelerationUpperLimit(const Robot& robot, 
                                                         const Eigen::VectorXd& amax)
  : ConstraintComponentBase(),
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


KinematicsLevel JointAccelerationUpperLimit::kinematicsLevel() const {
  return KinematicsLevel::AccelerationLevel;
}


bool JointAccelerationUpperLimit::isFeasible(Robot& robot, 
                                             const ContactStatus& contact_status,
                                             const GridInfo& grid_info,
                                             const SplitSolution& s,
                                             ConstraintComponentData& data) const {
  for (int i=0; i<dimc_; ++i) {
    if (s.a.tail(dimc_).coeff(i) > amax_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointAccelerationUpperLimit::setSlack(Robot& robot, 
                                           const ContactStatus& contact_status,
                                           const GridInfo& grid_info,
                                           const SplitSolution& s,
                                           ConstraintComponentData& data) const {
  data.slack = amax_ - s.a.tail(dimc_);
}


void JointAccelerationUpperLimit::evalConstraint(Robot& robot, 
                                                 const ContactStatus& contact_status,
                                                 const GridInfo& grid_info,
                                                 const SplitSolution& s,
                                                 ConstraintComponentData& data) const {
  data.residual = s.a.tail(dimc_) - amax_ + data.slack;
  computeComplementarySlackness(data);
  data.log_barrier = logBarrier(data.slack);
}


void JointAccelerationUpperLimit::evalDerivatives(
    Robot& robot, const ContactStatus& contact_status,
    const GridInfo& grid_info, const SplitSolution& s,
    ConstraintComponentData& data, SplitKKTResidual& kkt_residual) const {
  kkt_residual.la.tail(dimc_).noalias() += data.dual;
}


void JointAccelerationUpperLimit::condenseSlackAndDual(
    const ContactStatus& contact_status, const GridInfo& grid_info, 
    ConstraintComponentData& data, SplitKKTMatrix& kkt_matrix,
    SplitKKTResidual& kkt_residual) const {
  kkt_matrix.Qaa.diagonal().tail(dimc_).array()
      += data.dual.array() / data.slack.array();
  computeCondensingCoeffcient(data);
  kkt_residual.la.tail(dimc_).noalias() += data.cond;
}


void JointAccelerationUpperLimit::expandSlackAndDual(
    const ContactStatus& contact_status, const GridInfo& grid_info,
    const SplitDirection& d, ConstraintComponentData& data) const {
  data.dslack = - d.da().tail(dimc_) - data.residual;
  computeDualDirection(data);
}


int JointAccelerationUpperLimit::dimc() const {
  return dimc_;
}

} // namespace robotoc