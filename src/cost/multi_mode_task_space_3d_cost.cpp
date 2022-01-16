#include "robotoc/cost/multi_mode_task_space_3d_cost.hpp"


namespace robotoc {

MultiModeTaskSpace3DCost::MultiModeTaskSpace3DCost(const Robot& robot, 
                                                   const int frame_id)
  : CostFunctionComponentBase(),
    frame_id_(frame_id),
    x3d_ref_(),
    x3d_weight_(),
    x3df_weight_(Eigen::Vector3d::Zero()),
    x3di_weight_(Eigen::Vector3d::Zero()) {
  x3d_ref_[0] = Eigen::Vector3d::Zero();
  x3d_weight_[0] = Eigen::Vector3d::Zero();
}


MultiModeTaskSpace3DCost::MultiModeTaskSpace3DCost()
  : CostFunctionComponentBase(),
    frame_id_(0),
    x3d_ref_(),
    x3d_weight_(),
    x3df_weight_(),
    x3di_weight_() {
}


MultiModeTaskSpace3DCost::~MultiModeTaskSpace3DCost() {
}


void MultiModeTaskSpace3DCost::set_x3d_ref(const Eigen::Vector3d& x3d_ref,
                                           const int contact_mode_id) {
  x3d_ref_[contact_mode_id] = x3d_ref;
}


void MultiModeTaskSpace3DCost::set_x3d_ref(
    const Eigen::Vector3d& x3d_ref, const std::vector<int>& contact_mode_ids) {
  for (const int e : contact_mode_ids) {
    set_x3d_ref(x3d_ref, e);
  }
}


void MultiModeTaskSpace3DCost::set_x3d_weight(const Eigen::Vector3d& x3d_weight,
                                              const int contact_mode_id) {
  x3d_weight_[contact_mode_id] = x3d_weight;
}


void MultiModeTaskSpace3DCost::set_x3d_weight(
    const Eigen::Vector3d& x3d_weight, 
    const std::vector<int>& contact_mode_ids) {
  for (const int e : contact_mode_ids) {
    set_x3d_weight(x3d_weight, e);
  }
}


void MultiModeTaskSpace3DCost::set_x3df_weight(const Eigen::Vector3d& x3df_weight) {
  x3df_weight_ = x3df_weight;
}


void MultiModeTaskSpace3DCost::set_x3di_weight(const Eigen::Vector3d& x3di_weight) {
  x3di_weight_ = x3di_weight;
}


bool MultiModeTaskSpace3DCost::useKinematics() const {
  return true;
}


double MultiModeTaskSpace3DCost::evalStageCost(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const int time_stage, const double t, const double dt, 
    const SplitSolution& s) const {
  const int contact_mode_id = contact_status.contactModeId();
  const auto& x3d_ref = x3d_ref_.at(contact_mode_id);
  const auto& x3d_weight = x3d_weight_.at(contact_mode_id);
  double l = 0;
  data.diff_3d = robot.framePosition(frame_id_) - x3d_ref;
  l += (x3d_weight.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
  return 0.5 * dt * l;
}


void MultiModeTaskSpace3DCost::evalStageCostDerivatives(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const int time_stage, const double t, const double dt, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  const int contact_mode_id = contact_status.contactModeId();
  const auto& x3d_weight = x3d_weight_.at(contact_mode_id);
  data.J_6d.setZero();
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.J_3d.noalias() 
      = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
  kkt_residual.lq().noalias() 
      += dt * data.J_3d.transpose() * x3d_weight.asDiagonal() * data.diff_3d;
}


void MultiModeTaskSpace3DCost::evalStageCostHessian(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const int time_stage, const double t, const double dt, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  const int contact_mode_id = contact_status.contactModeId();
  const auto& x3d_weight = x3d_weight_.at(contact_mode_id);
  kkt_matrix.Qqq().noalias()
      += dt * data.J_3d.transpose() * x3d_weight.asDiagonal() * data.J_3d;
}


double MultiModeTaskSpace3DCost::evalTerminalCost(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s) const {
  const auto& x3d_ref = x3d_ref_.at(0);
  double l = 0;
  data.diff_3d = robot.framePosition(frame_id_) - x3d_ref;
  l += (x3df_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
  return 0.5 * l;
}


void MultiModeTaskSpace3DCost::evalTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  data.J_6d.setZero();
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.J_3d.noalias() 
      = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
  kkt_residual.lq().noalias() 
      += data.J_3d.transpose() * x3df_weight_.asDiagonal() * data.diff_3d;
}


void MultiModeTaskSpace3DCost::evalTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qqq().noalias()
      += data.J_3d.transpose() * x3df_weight_.asDiagonal() * data.J_3d;
}


double MultiModeTaskSpace3DCost::evalImpulseCost(
    Robot& robot,  const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const double t, const ImpulseSplitSolution& s) const {
  const int impulse_mode_id = impulse_status.impulseModeId();
  const auto& x3d_ref = x3d_ref_.at(impulse_mode_id);
  double l = 0;
  data.diff_3d = robot.framePosition(frame_id_) - x3d_ref;
  l += (x3di_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
  return 0.5 * l;
}


void MultiModeTaskSpace3DCost::evalImpulseCostDerivatives(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const double t, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  data.J_6d.setZero();
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.J_3d.noalias() 
      = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
  kkt_residual.lq().noalias() 
      += data.J_3d.transpose() * x3di_weight_.asDiagonal() * data.diff_3d;
}


void MultiModeTaskSpace3DCost::evalImpulseCostHessian(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const double t, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qqq().noalias()
      += data.J_3d.transpose() * x3di_weight_.asDiagonal() * data.J_3d;
}

} // namespace robotoc