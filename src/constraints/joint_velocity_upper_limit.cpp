#include "robotoc/constraints/joint_velocity_upper_limit.hpp"


namespace robotoc {

JointVelocityUpperLimit::JointVelocityUpperLimit(const Robot& robot)
  : ConstraintComponentBase(),
    dimc_(robot.jointVelocityLimit().size()),
    vmax_(robot.jointVelocityLimit()) {
}


JointVelocityUpperLimit::JointVelocityUpperLimit()
  : ConstraintComponentBase(),
    dimc_(0),
    vmax_() {
}


JointVelocityUpperLimit::~JointVelocityUpperLimit() {
}


KinematicsLevel JointVelocityUpperLimit::kinematicsLevel() const {
  return KinematicsLevel::VelocityLevel;
}


bool JointVelocityUpperLimit::isFeasible(Robot& robot, 
                                         const ContactStatus& contact_status, 
                                         const GridInfo& grid_info,
                                         const SplitSolution& s,
                                         ConstraintComponentData& data) const {
  for (int i=0; i<dimc_; ++i) {
    if (s.v.tail(dimc_).coeff(i) > vmax_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointVelocityUpperLimit::setSlack(Robot& robot, 
                                       const ContactStatus& contact_status, 
                                       const GridInfo& grid_info,
                                       const SplitSolution& s,
                                       ConstraintComponentData& data) const {
  data.slack = vmax_ - s.v.tail(dimc_);
}


void JointVelocityUpperLimit::evalConstraint(Robot& robot,
                                             const ContactStatus& contact_status,
                                             const GridInfo& grid_info,
                                             const SplitSolution& s,
                                             ConstraintComponentData& data) const {
  data.residual = s.v.tail(dimc_) - vmax_ + data.slack;
  computeComplementarySlackness(data);
  data.log_barrier = logBarrier(data.slack);
}


void JointVelocityUpperLimit::evalDerivatives(
    Robot& robot, const ContactStatus& contact_status,
    const GridInfo& grid_info, const SplitSolution& s,
    ConstraintComponentData& data, SplitKKTResidual& kkt_residual) const {
  kkt_residual.lv().tail(dimc_).noalias() += data.dual;
}


void JointVelocityUpperLimit::condenseSlackAndDual(
    const ContactStatus& contact_status, const GridInfo& grid_info, 
    ConstraintComponentData& data, SplitKKTMatrix& kkt_matrix,
    SplitKKTResidual& kkt_residual) const {
  kkt_matrix.Qvv().diagonal().tail(dimc_).array()
      += data.dual.array() / data.slack.array();
  computeCondensingCoeffcient(data);
  kkt_residual.lv().tail(dimc_).noalias() += data.cond;
}


void JointVelocityUpperLimit::expandSlackAndDual(
    const ContactStatus& contact_status, const GridInfo& grid_info,
    const SplitDirection& d, ConstraintComponentData& data) const {
  data.dslack = - d.dv().tail(dimc_) - data.residual;
  computeDualDirection(data);
}


int JointVelocityUpperLimit::dimc() const {
  return dimc_;
}

} // namespace robotoc