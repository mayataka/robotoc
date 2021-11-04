#include "robotoc/cost/time_varying_com_cost.hpp"


namespace robotoc {

TimeVaryingCoMCost::TimeVaryingCoMCost(
    const Robot& robot, 
    const std::shared_ptr<TimeVaryingCoMRefBase>& ref) 
  : CostFunctionComponentBase(),
    ref_(ref),
    q_weight_(Eigen::Vector3d::Zero()),
    qf_weight_(Eigen::Vector3d::Zero()),
    qi_weight_(Eigen::Vector3d::Zero()) {
}


TimeVaryingCoMCost::TimeVaryingCoMCost()
  : CostFunctionComponentBase(),
    ref_(),
    q_weight_(),
    qf_weight_(),
    qi_weight_() {
}


TimeVaryingCoMCost::~TimeVaryingCoMCost() {
}


void TimeVaryingCoMCost::set_ref(
    const std::shared_ptr<TimeVaryingCoMRefBase>& ref) {
  ref_ = ref;
}


void TimeVaryingCoMCost::set_q_weight(const Eigen::Vector3d& q_weight) {
  q_weight_ = q_weight;
}


void TimeVaryingCoMCost::set_qf_weight(const Eigen::Vector3d& qf_weight) {
  qf_weight_ = qf_weight;
}


void TimeVaryingCoMCost::set_qi_weight(const Eigen::Vector3d& qi_weight) {
  qi_weight_ = qi_weight;
}


bool TimeVaryingCoMCost::useKinematics() const {
  return true;
}


double TimeVaryingCoMCost::evalStageCost(Robot& robot, 
                                         const ContactStatus& contact_status, 
                                         CostFunctionData& data, 
                                         const double t, const double dt, 
                                         const SplitSolution& s) const {
  if (ref_->isActive(t)) {
    double l = 0;
    ref_->update_CoM_ref(t, data.q_3d_ref);
    data.diff_3d = robot.CoM() - data.q_3d_ref;
    l += (q_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
    return 0.5 * dt * l;
  }
  else {
    return 0;
  }
}


void TimeVaryingCoMCost::evalStageCostDerivatives(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const double t, const double dt, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) const {
  if (ref_->isActive(t)) {
    data.J_3d.setZero();
    robot.getCoMJacobian(data.J_3d);
    kkt_residual.lq().noalias() 
        += dt * data.J_3d.transpose() * q_weight_.asDiagonal() * data.diff_3d;
  }
}


void TimeVaryingCoMCost::evalStageCostHessian(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const double t, const double dt, const SplitSolution& s, 
    SplitKKTMatrix& kkt_matrix) const {
  if (ref_->isActive(t)) {
    kkt_matrix.Qqq().noalias()
        += dt * data.J_3d.transpose() * q_weight_.asDiagonal() * data.J_3d;
  }
}


double TimeVaryingCoMCost::evalTerminalCost(Robot& robot, 
                                            CostFunctionData& data, 
                                            const double t, 
                                            const SplitSolution& s) const {
  if (ref_->isActive(t)) {
    double l = 0;
    ref_->update_CoM_ref(t, data.q_3d_ref);
    data.diff_3d = robot.CoM() - data.q_3d_ref;
    l += (qf_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
    return 0.5 * l;
  }
  else {
    return 0;
  }
}


void TimeVaryingCoMCost::evalTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  if (ref_->isActive(t)) {
    data.J_3d.setZero();
    robot.getCoMJacobian(data.J_3d);
    kkt_residual.lq().noalias() 
        += data.J_3d.transpose() * qf_weight_.asDiagonal() * data.diff_3d;
  }
}


void TimeVaryingCoMCost::evalTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  if (ref_->isActive(t)) {
    kkt_matrix.Qqq().noalias()
        += data.J_3d.transpose() * qf_weight_.asDiagonal() * data.J_3d;
  }
}


double TimeVaryingCoMCost::evalImpulseCost(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const double t, const ImpulseSplitSolution& s) const {
  if (ref_->isActive(t)) {
    double l = 0;
    ref_->update_CoM_ref(t, data.q_3d_ref);
    data.diff_3d = robot.CoM() - data.q_3d_ref;
    l += (qi_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
    return 0.5 * l;
  }
  else {
    return 0;
  }
}


void TimeVaryingCoMCost::evalImpulseCostDerivatives(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const double t, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  if (ref_->isActive(t)) {
    data.J_3d.setZero();
    robot.getCoMJacobian(data.J_3d);
    kkt_residual.lq().noalias() 
        += data.J_3d.transpose() * qi_weight_.asDiagonal() * data.diff_3d;
  }
}


void TimeVaryingCoMCost::evalImpulseCostHessian(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const double t, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTMatrix& kkt_matrix) const {
  if (ref_->isActive(t)) {
    kkt_matrix.Qqq().noalias()
        += data.J_3d.transpose() * qi_weight_.asDiagonal() * data.J_3d;
  }
}

} // namespace robotoc