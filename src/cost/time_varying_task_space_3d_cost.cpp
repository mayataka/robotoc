#include "robotoc/cost/time_varying_task_space_3d_cost.hpp"


namespace robotoc {

TimeVaryingTaskSpace3DCost::TimeVaryingTaskSpace3DCost(
    const Robot& robot, const int frame_id, 
    const std::shared_ptr<TimeVaryingTaskSpace3DRefBase>& x3d_ref) 
  : CostFunctionComponentBase(),
    frame_id_(frame_id),
    x3d_ref_(x3d_ref),
    x3d_weight_(Eigen::Vector3d::Zero()),
    x3df_weight_(Eigen::Vector3d::Zero()),
    x3di_weight_(Eigen::Vector3d::Zero()) {
}


TimeVaryingTaskSpace3DCost::TimeVaryingTaskSpace3DCost()
  : CostFunctionComponentBase(),
    frame_id_(),
    x3d_ref_(),
    x3d_weight_(),
    x3df_weight_(),
    x3di_weight_() {
}


TimeVaryingTaskSpace3DCost::~TimeVaryingTaskSpace3DCost() {
}


void TimeVaryingTaskSpace3DCost::set_x3d_ref(
    const std::shared_ptr<TimeVaryingTaskSpace3DRefBase>& x3d_ref) {
  x3d_ref_ = x3d_ref;
}


void TimeVaryingTaskSpace3DCost::set_x3d_weight(
    const Eigen::Vector3d& x3d_weight) {
  x3d_weight_ = x3d_weight;
}


void TimeVaryingTaskSpace3DCost::set_x3df_weight(
    const Eigen::Vector3d& x3df_weight) {
  x3df_weight_ = x3df_weight;
}


void TimeVaryingTaskSpace3DCost::set_x3di_weight(
    const Eigen::Vector3d& x3di_weight) {
  x3di_weight_ = x3di_weight;
}


bool TimeVaryingTaskSpace3DCost::useKinematics() const {
  return true;
}


double TimeVaryingTaskSpace3DCost::evalStageCost(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const GridInfo& grid_info, const SplitSolution& s) const {
  if (x3d_ref_->isActive(grid_info)) {
    double l = 0;
    x3d_ref_->update_x3d_ref(grid_info, data.x3d_ref);
    data.diff_3d = robot.framePosition(frame_id_) - data.x3d_ref;
    l += (x3d_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
    return 0.5 * grid_info.dt * l;
  }
  else {
    return 0;
  }
}


void TimeVaryingTaskSpace3DCost::evalStageCostDerivatives(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const GridInfo& grid_info, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) const {
  if (x3d_ref_->isActive(grid_info)) {
    data.J_6d.setZero();
    robot.getFrameJacobian(frame_id_, data.J_6d);
    data.J_3d.noalias() 
        = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
    kkt_residual.lq().noalias() 
        += grid_info.dt * data.J_3d.transpose() * x3d_weight_.asDiagonal() * data.diff_3d;
  }
}


void TimeVaryingTaskSpace3DCost::evalStageCostHessian(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const GridInfo& grid_info, const SplitSolution& s, 
    SplitKKTMatrix& kkt_matrix) const {
  if (x3d_ref_->isActive(grid_info)) {
    kkt_matrix.Qqq().noalias()
        += grid_info.dt * data.J_3d.transpose() * x3d_weight_.asDiagonal() * data.J_3d;
  }
}


double TimeVaryingTaskSpace3DCost::evalTerminalCost(
    Robot& robot, CostFunctionData& data, const GridInfo& grid_info, 
    const SplitSolution& s) const {
  if (x3d_ref_->isActive(grid_info)) {
    double l = 0;
    x3d_ref_->update_x3d_ref(grid_info, data.x3d_ref);
    data.diff_3d = robot.framePosition(frame_id_) - data.x3d_ref;
    l += (x3df_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
    return 0.5 * l;
  }
  else {
    return 0;
  }
}


void TimeVaryingTaskSpace3DCost::evalTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const GridInfo& grid_info, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  if (x3d_ref_->isActive(grid_info)) {
    data.J_6d.setZero();
    robot.getFrameJacobian(frame_id_, data.J_6d);
    data.J_3d.noalias() 
        = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
    kkt_residual.lq().noalias() 
        += data.J_3d.transpose() * x3df_weight_.asDiagonal() * data.diff_3d;
  }
}


void TimeVaryingTaskSpace3DCost::evalTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const GridInfo& grid_info, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  if (x3d_ref_->isActive(grid_info)) {
    kkt_matrix.Qqq().noalias()
        += data.J_3d.transpose() * x3df_weight_.asDiagonal() * data.J_3d;
  }
}


double TimeVaryingTaskSpace3DCost::evalImpulseCost(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s) const {
  if (x3d_ref_->isActive(grid_info)) {
    double l = 0;
    x3d_ref_->update_x3d_ref(grid_info, data.x3d_ref);
    data.diff_3d = robot.framePosition(frame_id_) - data.x3d_ref;
    l += (x3di_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
    return 0.5 * l;
  }
  else {
    return 0;
  }
}


void TimeVaryingTaskSpace3DCost::evalImpulseCostDerivatives(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  if (x3d_ref_->isActive(grid_info)) {
    data.J_6d.setZero();
    robot.getFrameJacobian(frame_id_, data.J_6d);
    data.J_3d.noalias() 
        = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
    kkt_residual.lq().noalias() 
        += data.J_3d.transpose() * x3di_weight_.asDiagonal() * data.diff_3d;
  }
}


void TimeVaryingTaskSpace3DCost::evalImpulseCostHessian(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTMatrix& kkt_matrix) const {
  if (x3d_ref_->isActive(grid_info)) {
    kkt_matrix.Qqq().noalias()
        += data.J_3d.transpose() * x3di_weight_.asDiagonal() * data.J_3d;
  }
}

} // namespace robotoc