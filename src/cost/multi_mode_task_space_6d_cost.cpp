#include "robotoc/cost/multi_mode_task_space_6d_cost.hpp"


namespace robotoc {

MultiModeTaskSpace6DCost::MultiModeTaskSpace6DCost(const Robot& robot, 
                                                   const int frame_id)
  : CostFunctionComponentBase(),
    frame_id_(frame_id),
    x6d_ref_(),
    x6d_ref_inv_(),
    x6d_weight_(), 
    x6df_weight_(Vector6d::Zero()), 
    x6di_weight_(Vector6d::Zero()) {
  x6d_ref_[0] = SE3(Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero());
  x6d_ref_inv_[0] = x6d_ref_[0].inverse();
  x6d_weight_[0] = Vector6d::Zero();
}


MultiModeTaskSpace6DCost::MultiModeTaskSpace6DCost()
  : CostFunctionComponentBase(),
    frame_id_(0),
    x6d_ref_(),
    x6d_ref_inv_(),
    x6d_weight_(), 
    x6df_weight_(Vector6d::Zero()), 
    x6di_weight_(Vector6d::Zero()) {
}


MultiModeTaskSpace6DCost::~MultiModeTaskSpace6DCost() {
}


void MultiModeTaskSpace6DCost::set_x6d_ref(const Eigen::Vector3d& trans_ref, 
                                           const Eigen::Matrix3d& rot_ref,
                                           const int contact_mode_id) {
  x6d_ref_[contact_mode_id] = SE3(rot_ref, trans_ref);
  x6d_ref_inv_[contact_mode_id] = x6d_ref_[contact_mode_id].inverse();
}


void MultiModeTaskSpace6DCost::set_x6d_ref(const Eigen::Vector3d& trans_ref, 
                                           const Eigen::Matrix3d& rot_ref,
                                           const std::vector<int>& contact_mode_ids) {
  for (const int e : contact_mode_ids) {
    set_x6d_ref(trans_ref, rot_ref, e);
  }
}


void MultiModeTaskSpace6DCost::set_x6d_weight(const Eigen::Vector3d& trans_weight, 
                                              const Eigen::Vector3d& rot_weight,
                                              const int contact_mode_id) {
  x6d_weight_[contact_mode_id].template head<3>() = rot_weight;
  x6d_weight_[contact_mode_id].template tail<3>() = trans_weight;
}


void MultiModeTaskSpace6DCost::set_x6d_weight(const Eigen::Vector3d& trans_weight, 
                                              const Eigen::Vector3d& rot_weight,
                                              const std::vector<int>& contact_mode_ids) {
  for (const int e : contact_mode_ids) {
    set_x6d_weight(trans_weight, rot_weight, e);
  }
}


void MultiModeTaskSpace6DCost::set_x6df_weight(const Eigen::Vector3d& trans_weight, 
                                               const Eigen::Vector3d& rot_weight) {
  x6df_weight_.template head<3>() = rot_weight;
  x6df_weight_.template tail<3>() = trans_weight;
}


void MultiModeTaskSpace6DCost::set_x6di_weight(const Eigen::Vector3d& trans_weight, 
                                               const Eigen::Vector3d& rot_weight) {
  x6di_weight_.template head<3>() = rot_weight;
  x6di_weight_.template tail<3>() = trans_weight;
}
 

bool MultiModeTaskSpace6DCost::useKinematics() const {
  return true;
}


double MultiModeTaskSpace6DCost::evalStageCost(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const GridInfo& grid_info, const SplitSolution& s) const {
  const int contact_mode_id = contact_status.contactModeId();
  const auto& x6d_ref_inv = x6d_ref_inv_.at(contact_mode_id);
  const auto& x6d_weight = x6d_weight_.at(contact_mode_id);
  double l = 0;
  data.diff_x6d = x6d_ref_inv * robot.framePlacement(frame_id_);
  data.diff_6d = Log6Map(data.diff_x6d);
  l += (x6d_weight.array()*data.diff_6d.array()*data.diff_6d.array()).sum();
  return 0.5 * grid_info.dt * l;
}


void MultiModeTaskSpace6DCost::evalStageCostDerivatives(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const GridInfo& grid_info, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) const {
  const int contact_mode_id = contact_status.contactModeId();
  const auto& x6d_weight = x6d_weight_.at(contact_mode_id);
  data.J_66.setZero();
  computeJLog6Map(data.diff_x6d, data.J_66);
  data.J_6d.setZero();
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.JJ_6d.noalias() = data.J_66 * data.J_6d;
  kkt_residual.lq().noalias() 
      += grid_info.dt * data.JJ_6d.transpose() * x6d_weight.asDiagonal() * data.diff_6d;
}


void MultiModeTaskSpace6DCost::evalStageCostHessian(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const GridInfo& grid_info, const SplitSolution& s, 
    SplitKKTMatrix& kkt_matrix) const {
  const int contact_mode_id = contact_status.contactModeId();
  const auto& x6d_weight = x6d_weight_.at(contact_mode_id);
  kkt_matrix.Qqq().noalias()
      += grid_info.dt * data.JJ_6d.transpose() * x6d_weight.asDiagonal() * data.JJ_6d;
}


double MultiModeTaskSpace6DCost::evalTerminalCost(
    Robot& robot, CostFunctionData& data, const GridInfo& grid_info, 
    const SplitSolution& s) const {
  const auto& x6d_ref_inv = x6d_ref_inv_.at(0);
  double l = 0;
  data.diff_x6d = x6d_ref_inv * robot.framePlacement(frame_id_);
  data.diff_6d = Log6Map(data.diff_x6d);
  l += (x6df_weight_.array()*data.diff_6d.array()*data.diff_6d.array()).sum();
  return 0.5 * l;
}


void MultiModeTaskSpace6DCost::evalTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const GridInfo& grid_info, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  data.J_66.setZero();
  computeJLog6Map(data.diff_x6d, data.J_66);
  data.J_6d.setZero();
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.JJ_6d.noalias() = data.J_66 * data.J_6d;
  kkt_residual.lq().noalias() 
      += data.JJ_6d.transpose() * x6df_weight_.asDiagonal() * data.diff_6d;
}


void MultiModeTaskSpace6DCost::evalTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const GridInfo& grid_info, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qqq().noalias()
      += data.JJ_6d.transpose() * x6df_weight_.asDiagonal() * data.JJ_6d;
}


double MultiModeTaskSpace6DCost::evalImpulseCost(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s) const {
  const int impulse_mode_id = impulse_status.impulseModeId();
  const auto& x6d_ref_inv = x6d_ref_inv_.at(impulse_mode_id);
  double l = 0;
  data.diff_x6d = x6d_ref_inv * robot.framePlacement(frame_id_);
  data.diff_6d = Log6Map(data.diff_x6d);
  l += (x6di_weight_.array()*data.diff_6d.array()*data.diff_6d.array()).sum();
  return 0.5 * l;
}


void MultiModeTaskSpace6DCost::evalImpulseCostDerivatives(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  data.J_66.setZero();
  computeJLog6Map(data.diff_x6d, data.J_66);
  data.J_6d.setZero();
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.JJ_6d.noalias() = data.J_66 * data.J_6d;
  kkt_residual.lq().noalias() 
      += data.JJ_6d.transpose() * x6di_weight_.asDiagonal() * data.diff_6d;
}


void MultiModeTaskSpace6DCost::evalImpulseCostHessian(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qqq().noalias()
      += data.JJ_6d.transpose() * x6di_weight_.asDiagonal() * data.JJ_6d;
}

} // namespace robotoc