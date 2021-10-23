#include "robotoc/cost/task_space_6d_cost.hpp"


namespace robotoc {

TaskSpace6DCost::TaskSpace6DCost(const Robot& robot, const int frame_id)
  : CostFunctionComponentBase(),
    frame_id_(frame_id),
    SE3_ref_(SE3(Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero())),
    SE3_ref_inv_(SE3_ref_.inverse()),
    q_6d_weight_(Eigen::VectorXd::Zero(6)), 
    qf_6d_weight_(Eigen::VectorXd::Zero(6)), 
    qi_6d_weight_(Eigen::VectorXd::Zero(6)) {
}


TaskSpace6DCost::TaskSpace6DCost()
  : CostFunctionComponentBase(),
    frame_id_(0),
    SE3_ref_(SE3(Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero())),
    SE3_ref_inv_(SE3_ref_.inverse()),
    q_6d_weight_(Eigen::VectorXd::Zero(6)), 
    qf_6d_weight_(Eigen::VectorXd::Zero(6)), 
    qi_6d_weight_(Eigen::VectorXd::Zero(6)) {
}


TaskSpace6DCost::~TaskSpace6DCost() {
}


void TaskSpace6DCost::set_q_6d_ref(const Eigen::Vector3d& position_ref, 
                                   const Eigen::Matrix3d& rotation_mat_ref) {
  SE3_ref_ = SE3(rotation_mat_ref, position_ref);
  SE3_ref_inv_ = SE3_ref_.inverse();
}


void TaskSpace6DCost::set_q_weight(const Eigen::Vector3d& position_weight, 
                                   const Eigen::Vector3d& rotation_weight) {
  q_6d_weight_.template head<3>() = rotation_weight;
  q_6d_weight_.template tail<3>() = position_weight;
}


void TaskSpace6DCost::set_qf_weight(const Eigen::Vector3d& position_weight, 
                                    const Eigen::Vector3d& rotation_weight) {
  qf_6d_weight_.template head<3>() = rotation_weight;
  qf_6d_weight_.template tail<3>() = position_weight;
}


void TaskSpace6DCost::set_qi_weight(const Eigen::Vector3d& position_weight, 
                                    const Eigen::Vector3d& rotation_weight) {
  qi_6d_weight_.template head<3>() = rotation_weight;
  qi_6d_weight_.template tail<3>() = position_weight;
}
 

bool TaskSpace6DCost::useKinematics() const {
  return true;
}


double TaskSpace6DCost::computeStageCost(Robot& robot, CostFunctionData& data, 
                                         const double t, const double dt, 
                                         const SplitSolution& s) const {
  double l = 0;
  data.diff_SE3 = SE3_ref_inv_ * robot.framePlacement(frame_id_);
  data.diff_6d = Log6Map(data.diff_SE3);
  l += (q_6d_weight_.array()*data.diff_6d.array()*data.diff_6d.array()).sum();
  return 0.5 * dt * l;
}


void TaskSpace6DCost::computeStageCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, const double dt, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  data.J_66.setZero();
  computeJLog6Map(data.diff_SE3, data.J_66);
  data.J_6d.setZero();
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.JJ_6d.noalias() = data.J_66 * data.J_6d;
  kkt_residual.lq().noalias() 
      += dt * data.JJ_6d.transpose() * q_6d_weight_.asDiagonal() 
                                       * data.diff_6d;
}


void TaskSpace6DCost::computeStageCostHessian(
    Robot& robot, CostFunctionData& data, const double t, const double dt, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qqq().noalias()
      += dt * data.JJ_6d.transpose() * q_6d_weight_.asDiagonal() * data.JJ_6d;
}


double TaskSpace6DCost::computeTerminalCost(Robot& robot, 
                                            CostFunctionData& data, 
                                            const double t, 
                                            const SplitSolution& s) const {
  double l = 0;
  data.diff_SE3 = SE3_ref_inv_ * robot.framePlacement(frame_id_);
  data.diff_6d = Log6Map(data.diff_SE3);
  l += (qf_6d_weight_.array()*data.diff_6d.array()*data.diff_6d.array()).sum();
  return 0.5 * l;
}


void TaskSpace6DCost::computeTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  data.J_66.setZero();
  computeJLog6Map(data.diff_SE3, data.J_66);
  data.J_6d.setZero();
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.JJ_6d.noalias() = data.J_66 * data.J_6d;
  kkt_residual.lq().noalias() 
      += data.JJ_6d.transpose() * qf_6d_weight_.asDiagonal() * data.diff_6d;
}


void TaskSpace6DCost::computeTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qqq().noalias()
      += data.JJ_6d.transpose() * qf_6d_weight_.asDiagonal() * data.JJ_6d;
}


double TaskSpace6DCost::computeImpulseCost(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s) const {
  double l = 0;
  data.diff_SE3 = SE3_ref_inv_ * robot.framePlacement(frame_id_);
  data.diff_6d = Log6Map(data.diff_SE3);
  l += (qi_6d_weight_.array()*data.diff_6d.array()*data.diff_6d.array()).sum();
  return 0.5 * l;
}


void TaskSpace6DCost::computeImpulseCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  data.J_66.setZero();
  computeJLog6Map(data.diff_SE3, data.J_66);
  data.J_6d.setZero();
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.JJ_6d.noalias() = data.J_66 * data.J_6d;
  kkt_residual.lq().noalias() 
      += data.JJ_6d.transpose() * qi_6d_weight_.asDiagonal() * data.diff_6d;
}


void TaskSpace6DCost::computeImpulseCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qqq().noalias()
      += data.JJ_6d.transpose() * qi_6d_weight_.asDiagonal() * data.JJ_6d;
}

} // namespace robotoc