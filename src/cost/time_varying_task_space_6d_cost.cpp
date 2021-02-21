#include "idocp/cost/time_varying_task_space_6d_cost.hpp"

#include <iostream>


namespace idocp {

TimeVaryingTaskSpace6DCost::TimeVaryingTaskSpace6DCost(
    const Robot& robot, const int frame_id, 
    const std::shared_ptr<TimeVaryingTaskSpace6DRefBase>& ref) 
  : CostFunctionComponentBase(),
    frame_id_(frame_id),
    ref_(ref),
    q_6d_weight_(Eigen::VectorXd::Zero(6)), 
    qf_6d_weight_(Eigen::VectorXd::Zero(6)), 
    qi_6d_weight_(Eigen::VectorXd::Zero(6)) {
}


TimeVaryingTaskSpace6DCost::TimeVaryingTaskSpace6DCost()
  : CostFunctionComponentBase(),
    frame_id_(0),
    ref_(),
    q_6d_weight_(Eigen::VectorXd::Zero(6)), 
    qf_6d_weight_(Eigen::VectorXd::Zero(6)), 
    qi_6d_weight_(Eigen::VectorXd::Zero(6)) {
}


TimeVaryingTaskSpace6DCost::~TimeVaryingTaskSpace6DCost() {
}


bool TimeVaryingTaskSpace6DCost::useKinematics() const {
  return true;
}


void TimeVaryingTaskSpace6DCost::set_ref(
    const std::shared_ptr<TimeVaryingTaskSpace6DRefBase>& ref) {
  ref_ = ref;
}


void TimeVaryingTaskSpace6DCost::set_q_6d_weight(
    const Eigen::Vector3d& position_weight, 
    const Eigen::Vector3d& rotation_weight) {
  q_6d_weight_.template head<3>() = rotation_weight;
  q_6d_weight_.template tail<3>() = position_weight;
}


void TimeVaryingTaskSpace6DCost::set_qf_6d_weight(
    const Eigen::Vector3d& position_weight, 
    const Eigen::Vector3d& rotation_weight) {
  qf_6d_weight_.template head<3>() = rotation_weight;
  qf_6d_weight_.template tail<3>() = position_weight;
}


void TimeVaryingTaskSpace6DCost::set_qi_6d_weight(
    const Eigen::Vector3d& position_weight, 
    const Eigen::Vector3d& rotation_weight) {
  qi_6d_weight_.template head<3>() = rotation_weight;
  qi_6d_weight_.template tail<3>() = position_weight;
}
 

double TimeVaryingTaskSpace6DCost::computeStageCost(
    Robot& robot, CostFunctionData& data, const double t, const double dt, 
    const SplitSolution& s) const {
  double l = 0;
  ref_->compute_q_6d_ref(t, data.SE3_ref);
  data.SE3_ref_inv = data.SE3_ref.inverse();
  data.diff_SE3 = data.SE3_ref_inv * robot.framePlacement(frame_id_);
  data.diff_6d = pinocchio::log6(data.diff_SE3).toVector();
  l += (q_6d_weight_.array()*data.diff_6d.array()*data.diff_6d.array()).sum();
  return 0.5 * dt * l;
}


double TimeVaryingTaskSpace6DCost::computeTerminalCost(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s) const {
  double l = 0;
  ref_->compute_q_6d_ref(t, data.SE3_ref);
  data.SE3_ref_inv = data.SE3_ref.inverse();
  data.diff_SE3 = data.SE3_ref_inv * robot.framePlacement(frame_id_);
  data.diff_6d = pinocchio::log6(data.diff_SE3).toVector();
  l += (qf_6d_weight_.array()*data.diff_6d.array()*data.diff_6d.array()).sum();
  return 0.5 * l;
}


double TimeVaryingTaskSpace6DCost::computeImpulseCost(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s) const {
  double l = 0;
  ref_->compute_q_6d_ref(t, data.SE3_ref);
  data.SE3_ref_inv = data.SE3_ref.inverse();
  data.diff_SE3 = data.SE3_ref_inv * robot.framePlacement(frame_id_);
  data.diff_6d = pinocchio::log6(data.diff_SE3).toVector();
  l += (qi_6d_weight_.array()*data.diff_6d.array()*data.diff_6d.array()).sum();
  return 0.5 * l;
}


void TimeVaryingTaskSpace6DCost::computeStageCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, const double dt, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  ref_->compute_q_6d_ref(t, data.SE3_ref);
  data.SE3_ref_inv = data.SE3_ref.inverse();
  data.diff_SE3 = data.SE3_ref_inv * robot.framePlacement(frame_id_);
  data.diff_6d = pinocchio::log6(data.diff_SE3).toVector();
  pinocchio::Jlog6(data.diff_SE3, data.J_66);
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.JJ_6d.noalias() = data.J_66 * data.J_6d;
  kkt_residual.lq().noalias() 
      += dt * data.JJ_6d.transpose() * q_6d_weight_.asDiagonal() 
                                       * data.diff_6d;
}


void TimeVaryingTaskSpace6DCost::computeTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  ref_->compute_q_6d_ref(t, data.SE3_ref);
  data.SE3_ref_inv = data.SE3_ref.inverse();
  data.diff_SE3 = data.SE3_ref_inv * robot.framePlacement(frame_id_);
  data.diff_6d = pinocchio::log6(data.diff_SE3).toVector();
  pinocchio::Jlog6(data.diff_SE3, data.J_66);
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.JJ_6d.noalias() = data.J_66 * data.J_6d;
  kkt_residual.lq().noalias() 
      += data.JJ_6d.transpose() * qf_6d_weight_.asDiagonal() * data.diff_6d;
}


void TimeVaryingTaskSpace6DCost::computeImpulseCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  ref_->compute_q_6d_ref(t, data.SE3_ref);
  data.SE3_ref_inv = data.SE3_ref.inverse();
  data.diff_SE3 = data.SE3_ref_inv * robot.framePlacement(frame_id_);
  data.diff_6d = pinocchio::log6(data.diff_SE3).toVector();
  pinocchio::Jlog6(data.diff_SE3, data.J_66);
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.JJ_6d.noalias() = data.J_66 * data.J_6d;
  kkt_residual.lq().noalias() 
      += data.JJ_6d.transpose() * qi_6d_weight_.asDiagonal() * data.diff_6d;
}


void TimeVaryingTaskSpace6DCost::computeStageCostHessian(
    Robot& robot, CostFunctionData& data, const double t, const double dt, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  ref_->compute_q_6d_ref(t, data.SE3_ref);
  data.SE3_ref_inv = data.SE3_ref.inverse();
  data.diff_SE3 = data.SE3_ref_inv * robot.framePlacement(frame_id_);
  pinocchio::Jlog6(data.diff_SE3, data.J_66);
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.JJ_6d.noalias() = data.J_66 * data.J_6d;
  kkt_matrix.Qqq().noalias()
      += dt * data.JJ_6d.transpose() * q_6d_weight_.asDiagonal() * data.JJ_6d;
}


void TimeVaryingTaskSpace6DCost::computeTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  ref_->compute_q_6d_ref(t, data.SE3_ref);
  data.SE3_ref_inv = data.SE3_ref.inverse();
  data.diff_SE3 = data.SE3_ref_inv * robot.framePlacement(frame_id_);
  pinocchio::Jlog6(data.diff_SE3, data.J_66);
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.JJ_6d.noalias() = data.J_66 * data.J_6d;
  kkt_matrix.Qqq().noalias()
      += data.JJ_6d.transpose() * qf_6d_weight_.asDiagonal() * data.JJ_6d;
}


void TimeVaryingTaskSpace6DCost::computeImpulseCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTMatrix& kkt_matrix) const {
  ref_->compute_q_6d_ref(t, data.SE3_ref);
  data.SE3_ref_inv = data.SE3_ref.inverse();
  data.diff_SE3 = data.SE3_ref_inv * robot.framePlacement(frame_id_);
  pinocchio::Jlog6(data.diff_SE3, data.J_66);
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.JJ_6d.noalias() = data.J_66 * data.J_6d;
  kkt_matrix.Qqq().noalias()
      += data.JJ_6d.transpose() * qi_6d_weight_.asDiagonal() * data.JJ_6d;
}

} // namespace idocp