#include "robotoc/constraints/joint_acceleration_lower_limit.hpp"


namespace robotoc {

JointAccelerationLowerLimit::JointAccelerationLowerLimit(const Robot& robot, 
                                                         const Eigen::VectorXd& amin)
  : ConstraintComponentBase(),
    dimc_(amin.size()),
    amin_(amin) {
}


JointAccelerationLowerLimit::JointAccelerationLowerLimit()
  : ConstraintComponentBase(),
    dimc_(0),
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
                                             const ContactStatus& contact_status,
                                             ConstraintComponentData& data, 
                                             const SplitSolution& s) const {
  for (int i=0; i<dimc_; ++i) {
    if (s.a.tail(dimc_).coeff(i) < amin_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointAccelerationLowerLimit::setSlack(Robot& robot, 
                                           const ContactStatus& contact_status,
                                           ConstraintComponentData& data, 
                                           const SplitSolution& s) const {
  data.slack = s.a.tail(dimc_) - amin_;
}


void JointAccelerationLowerLimit::evalConstraint(Robot& robot, 
                                                 const ContactStatus& contact_status,
                                                 ConstraintComponentData& data, 
                                                 const SplitSolution& s) const {
  data.residual = amin_ - s.a.tail(dimc_) + data.slack;
  computeComplementarySlackness(data);
  data.log_barrier = logBarrier(data.slack);
}


void JointAccelerationLowerLimit::evalDerivatives(
    Robot& robot, const ContactStatus& contact_status,
    ConstraintComponentData& data, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) const {
  kkt_residual.la.tail(dimc_).noalias() -= data.dual;
}


void JointAccelerationLowerLimit::condenseSlackAndDual(
    const ContactStatus& contact_status, ConstraintComponentData& data, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) const {
  kkt_matrix.Qaa.diagonal().tail(dimc_).array()
      += data.dual.array() / data.slack.array();
  computeCondensingCoeffcient(data);
  kkt_residual.la.tail(dimc_).noalias() -= data.cond;
}


void JointAccelerationLowerLimit::expandSlackAndDual(
    const ContactStatus& contact_status, ConstraintComponentData& data, 
    const SplitDirection& d) const {
  data.dslack = d.da().tail(dimc_) - data.residual;
  computeDualDirection(data);
}


int JointAccelerationLowerLimit::dimc() const {
  return dimc_;
}

} // namespace robotoc