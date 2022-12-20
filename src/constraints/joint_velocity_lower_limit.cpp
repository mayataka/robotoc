#include "robotoc/constraints/joint_velocity_lower_limit.hpp"


namespace robotoc {

JointVelocityLowerLimit::JointVelocityLowerLimit(const Robot& robot)
  : ConstraintComponentBase(),
    dimc_(robot.jointVelocityLimit().size()),
    vmin_(-robot.jointVelocityLimit()) {
}


JointVelocityLowerLimit::JointVelocityLowerLimit()
  : ConstraintComponentBase(),
    dimc_(0),
    vmin_() {
}


JointVelocityLowerLimit::~JointVelocityLowerLimit() {
}


KinematicsLevel JointVelocityLowerLimit::kinematicsLevel() const {
  return KinematicsLevel::VelocityLevel;
}


bool JointVelocityLowerLimit::isFeasible(Robot& robot, 
                                        const ContactStatus& contact_status, 
                                        const GridInfo& grid_info,
                                        const SplitSolution& s,
                                        ConstraintComponentData& data) const {
  for (int i=0; i<dimc_; ++i) {
    if (s.v.tail(dimc_).coeff(i) < vmin_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointVelocityLowerLimit::setSlack(Robot& robot, 
                                       const ContactStatus& contact_status, 
                                       const GridInfo& grid_info,
                                       const SplitSolution& s,
                                       ConstraintComponentData& data) const {
  data.slack = s.v.tail(dimc_) - vmin_;
}


void JointVelocityLowerLimit::evalConstraint(Robot& robot, 
                                             const ContactStatus& contact_status,
                                             const GridInfo& grid_info,
                                             const SplitSolution& s,
                                             ConstraintComponentData& data) const {
  data.residual = vmin_ - s.v.tail(dimc_) + data.slack;
  computeComplementarySlackness(data);
  data.log_barrier = logBarrier(data.slack);
}


void JointVelocityLowerLimit::evalDerivatives(
    Robot& robot, const ContactStatus& contact_status,
    const GridInfo& grid_info, const SplitSolution& s,
    ConstraintComponentData& data, SplitKKTResidual& kkt_residual) const {
  kkt_residual.lv().tail(dimc_).noalias() -= data.dual;
}


void JointVelocityLowerLimit::condenseSlackAndDual(
    const ContactStatus& contact_status, const GridInfo& grid_info, 
    ConstraintComponentData& data, SplitKKTMatrix& kkt_matrix,
    SplitKKTResidual& kkt_residual) const {
  kkt_matrix.Qvv().diagonal().tail(dimc_).array()
      += data.dual.array() / data.slack.array();
  computeCondensingCoeffcient(data);
  kkt_residual.lv().tail(dimc_).noalias() -= data.cond;
}


void JointVelocityLowerLimit::expandSlackAndDual(
    const ContactStatus& contact_status, const GridInfo& grid_info,
    const SplitDirection& d, ConstraintComponentData& data) const {
  data.dslack = d.dv().tail(dimc_) - data.residual;
  computeDualDirection(data);
}


int JointVelocityLowerLimit::dimc() const {
  return dimc_;
}

} // namespace robotoc