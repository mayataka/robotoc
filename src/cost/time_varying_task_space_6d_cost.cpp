#include "robotoc/cost/time_varying_task_space_6d_cost.hpp"


namespace robotoc {

TimeVaryingTaskSpace6DCost::TimeVaryingTaskSpace6DCost(
    const Robot& robot, const int frame_id, 
    const std::shared_ptr<TimeVaryingTaskSpace6DRefBase>& x6d_ref) 
  : CostFunctionComponentBase(),
    frame_id_(frame_id),
    x6d_ref_(x6d_ref),
    x6d_weight_(Eigen::VectorXd::Zero(6)), 
    x6df_weight_(Eigen::VectorXd::Zero(6)), 
    x6di_weight_(Eigen::VectorXd::Zero(6)) {
}


TimeVaryingTaskSpace6DCost::TimeVaryingTaskSpace6DCost()
  : CostFunctionComponentBase(),
    frame_id_(0),
    x6d_ref_(),
    x6d_weight_(Eigen::VectorXd::Zero(6)), 
    x6df_weight_(Eigen::VectorXd::Zero(6)), 
    x6di_weight_(Eigen::VectorXd::Zero(6)) {
}


TimeVaryingTaskSpace6DCost::~TimeVaryingTaskSpace6DCost() {
}


void TimeVaryingTaskSpace6DCost::set_x6d_ref(
    const std::shared_ptr<TimeVaryingTaskSpace6DRefBase>& x6d_ref) {
  x6d_ref_ = x6d_ref;
}


void TimeVaryingTaskSpace6DCost::set_x6d_weight(
    const Eigen::Vector3d& trans_weight, 
    const Eigen::Vector3d& rot_weight) {
  x6d_weight_.template head<3>() = rot_weight;
  x6d_weight_.template tail<3>() = trans_weight;
}


void TimeVaryingTaskSpace6DCost::set_x6df_weight(
    const Eigen::Vector3d& trans_weight, 
    const Eigen::Vector3d& rot_weight) {
  x6df_weight_.template head<3>() = rot_weight;
  x6df_weight_.template tail<3>() = trans_weight;
}


void TimeVaryingTaskSpace6DCost::set_x6di_weight(
    const Eigen::Vector3d& trans_weight, 
    const Eigen::Vector3d& rot_weight) {
  x6di_weight_.template head<3>() = rot_weight;
  x6di_weight_.template tail<3>() = trans_weight;
}
 

bool TimeVaryingTaskSpace6DCost::useKinematics() const {
  return true;
}


double TimeVaryingTaskSpace6DCost::evalStageCost(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const GridInfo& grid_info, const SplitSolution& s) const {
  if (x6d_ref_->isActive(grid_info.t)) {
    double l = 0;
    x6d_ref_->update_x6d_ref(grid_info.t, data.x6d_ref);
    data.x6d_ref_inv = data.x6d_ref.inverse();
    data.diff_x6d = data.x6d_ref_inv * robot.framePlacement(frame_id_);
    data.diff_6d = Log6Map(data.diff_x6d);
    l += (x6d_weight_.array()*data.diff_6d.array()*data.diff_6d.array()).sum();
    return 0.5 * grid_info.dt * l;
  }
  else {
    return 0;
  }
}


void TimeVaryingTaskSpace6DCost::evalStageCostDerivatives(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const GridInfo& grid_info, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) const {
  if (x6d_ref_->isActive(grid_info.t)) {
    data.J_66.setZero();
    computeJLog6Map(data.diff_x6d, data.J_66);
    data.J_6d.setZero();
    robot.getFrameJacobian(frame_id_, data.J_6d);
    data.JJ_6d.noalias() = data.J_66 * data.J_6d;
    kkt_residual.lq().noalias() 
        += grid_info.dt * data.JJ_6d.transpose() * x6d_weight_.asDiagonal() * data.diff_6d;
  }
}


void TimeVaryingTaskSpace6DCost::evalStageCostHessian(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const GridInfo& grid_info, const SplitSolution& s, 
    SplitKKTMatrix& kkt_matrix) const {
  if (x6d_ref_->isActive(grid_info.t)) {
    kkt_matrix.Qqq().noalias()
        += grid_info.dt * data.JJ_6d.transpose() * x6d_weight_.asDiagonal() * data.JJ_6d;
  }
}


double TimeVaryingTaskSpace6DCost::evalTerminalCost(
    Robot& robot, CostFunctionData& data, const GridInfo& grid_info, 
    const SplitSolution& s) const {
  if (x6d_ref_->isActive(grid_info.t)) {
    double l = 0;
    x6d_ref_->update_x6d_ref(grid_info.t, data.x6d_ref);
    data.x6d_ref_inv = data.x6d_ref.inverse();
    data.diff_x6d = data.x6d_ref_inv * robot.framePlacement(frame_id_);
    data.diff_6d = Log6Map(data.diff_x6d);
    l += (x6df_weight_.array()*data.diff_6d.array()*data.diff_6d.array()).sum();
    return 0.5 * l;
  }
  else {
    return 0;
  }
}


void TimeVaryingTaskSpace6DCost::evalTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const GridInfo& grid_info, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  if (x6d_ref_->isActive(grid_info.t)) {
    data.J_66.setZero();
    computeJLog6Map(data.diff_x6d, data.J_66);
    data.J_6d.setZero();
    robot.getFrameJacobian(frame_id_, data.J_6d);
    data.JJ_6d.noalias() = data.J_66 * data.J_6d;
    kkt_residual.lq().noalias() 
        += data.JJ_6d.transpose() * x6df_weight_.asDiagonal() * data.diff_6d;
  }
}


void TimeVaryingTaskSpace6DCost::evalTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const GridInfo& grid_info, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  if (x6d_ref_->isActive(grid_info.t)) {
    kkt_matrix.Qqq().noalias()
        += data.JJ_6d.transpose() * x6df_weight_.asDiagonal() * data.JJ_6d;
  }
}


double TimeVaryingTaskSpace6DCost::evalImpulseCost(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s) const {
  if (x6d_ref_->isActive(grid_info.t)) {
    double l = 0;
    x6d_ref_->update_x6d_ref(grid_info.t, data.x6d_ref);
    data.x6d_ref_inv = data.x6d_ref.inverse();
    data.diff_x6d = data.x6d_ref_inv * robot.framePlacement(frame_id_);
    data.diff_6d = Log6Map(data.diff_x6d);
    l += (x6di_weight_.array()*data.diff_6d.array()*data.diff_6d.array()).sum();
    return 0.5 * l;
  }
  else {
    return 0;
  }
}


void TimeVaryingTaskSpace6DCost::evalImpulseCostDerivatives(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  if (x6d_ref_->isActive(grid_info.t)) {
    data.J_66.setZero();
    computeJLog6Map(data.diff_x6d, data.J_66);
    data.J_6d.setZero();
    robot.getFrameJacobian(frame_id_, data.J_6d);
    data.JJ_6d.noalias() = data.J_66 * data.J_6d;
    kkt_residual.lq().noalias() 
        += data.JJ_6d.transpose() * x6di_weight_.asDiagonal() * data.diff_6d;
  }
}


void TimeVaryingTaskSpace6DCost::evalImpulseCostHessian(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTMatrix& kkt_matrix) const {
  if (x6d_ref_->isActive(grid_info.t)) {
    kkt_matrix.Qqq().noalias()
        += data.JJ_6d.transpose() * x6di_weight_.asDiagonal() * data.JJ_6d;
  }
}

} // namespace robotoc