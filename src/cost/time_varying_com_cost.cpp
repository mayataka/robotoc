#include "robotoc/cost/time_varying_com_cost.hpp"


namespace robotoc {

TimeVaryingCoMCost::TimeVaryingCoMCost(
    const Robot& robot, 
    const std::shared_ptr<TimeVaryingCoMRefBase>& com_ref) 
  : CostFunctionComponentBase(),
    com_ref_(com_ref),
    com_weight_(Eigen::Vector3d::Zero()),
    comf_weight_(Eigen::Vector3d::Zero()),
    comi_weight_(Eigen::Vector3d::Zero()) {
}


TimeVaryingCoMCost::TimeVaryingCoMCost()
  : CostFunctionComponentBase(),
    com_ref_(),
    com_weight_(),
    comf_weight_(),
    comi_weight_() {
}


TimeVaryingCoMCost::~TimeVaryingCoMCost() {
}


void TimeVaryingCoMCost::set_com_ref(
    const std::shared_ptr<TimeVaryingCoMRefBase>& com_ref) {
  com_ref_ = com_ref;
}


void TimeVaryingCoMCost::set_com_weight(const Eigen::Vector3d& com_weight) {
  com_weight_ = com_weight;
}


void TimeVaryingCoMCost::set_comf_weight(const Eigen::Vector3d& comf_weight) {
  comf_weight_ = comf_weight;
}


void TimeVaryingCoMCost::set_comi_weight(const Eigen::Vector3d& comi_weight) {
  comi_weight_ = comi_weight;
}


bool TimeVaryingCoMCost::useKinematics() const {
  return true;
}


double TimeVaryingCoMCost::evalStageCost(Robot& robot, 
                                         const ContactStatus& contact_status, 
                                         CostFunctionData& data, 
                                         const int time_stage,
                                         const double t, const double dt, 
                                         const SplitSolution& s) const {
  if (com_ref_->isActive(t)) {
    double l = 0;
    com_ref_->update_com_ref(t, data.x3d_ref);
    data.diff_3d = robot.CoM() - data.x3d_ref;
    l += (com_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
    return 0.5 * dt * l;
  }
  else {
    return 0;
  }
}


void TimeVaryingCoMCost::evalStageCostDerivatives(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const int time_stage, const double t, const double dt, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  if (com_ref_->isActive(t)) {
    data.J_3d.setZero();
    robot.getCoMJacobian(data.J_3d);
    kkt_residual.lq().noalias() 
        += dt * data.J_3d.transpose() * com_weight_.asDiagonal() * data.diff_3d;
  }
}


void TimeVaryingCoMCost::evalStageCostHessian(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const int time_stage, const double t, const double dt, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  if (com_ref_->isActive(t)) {
    kkt_matrix.Qqq().noalias()
        += dt * data.J_3d.transpose() * com_weight_.asDiagonal() * data.J_3d;
  }
}


double TimeVaryingCoMCost::evalTerminalCost(Robot& robot, 
                                            CostFunctionData& data, 
                                            const double t, 
                                            const SplitSolution& s) const {
  if (com_ref_->isActive(t)) {
    double l = 0;
    com_ref_->update_com_ref(t, data.x3d_ref);
    data.diff_3d = robot.CoM() - data.x3d_ref;
    l += (comf_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
    return 0.5 * l;
  }
  else {
    return 0;
  }
}


void TimeVaryingCoMCost::evalTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  if (com_ref_->isActive(t)) {
    data.J_3d.setZero();
    robot.getCoMJacobian(data.J_3d);
    kkt_residual.lq().noalias() 
        += data.J_3d.transpose() * comf_weight_.asDiagonal() * data.diff_3d;
  }
}


void TimeVaryingCoMCost::evalTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  if (com_ref_->isActive(t)) {
    kkt_matrix.Qqq().noalias()
        += data.J_3d.transpose() * comf_weight_.asDiagonal() * data.J_3d;
  }
}


double TimeVaryingCoMCost::evalImpulseCost(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const double t, const ImpulseSplitSolution& s) const {
  if (com_ref_->isActive(t)) {
    double l = 0;
    com_ref_->update_com_ref(t, data.x3d_ref);
    data.diff_3d = robot.CoM() - data.x3d_ref;
    l += (comi_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
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
  if (com_ref_->isActive(t)) {
    data.J_3d.setZero();
    robot.getCoMJacobian(data.J_3d);
    kkt_residual.lq().noalias() 
        += data.J_3d.transpose() * comi_weight_.asDiagonal() * data.diff_3d;
  }
}


void TimeVaryingCoMCost::evalImpulseCostHessian(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const double t, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTMatrix& kkt_matrix) const {
  if (com_ref_->isActive(t)) {
    kkt_matrix.Qqq().noalias()
        += data.J_3d.transpose() * comi_weight_.asDiagonal() * data.J_3d;
  }
}

} // namespace robotoc