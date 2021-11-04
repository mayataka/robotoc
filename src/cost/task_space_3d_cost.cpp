#include "robotoc/cost/task_space_3d_cost.hpp"


namespace robotoc {

TaskSpace3DCost::TaskSpace3DCost(const Robot& robot, const int frame_id)
  : CostFunctionComponentBase(),
    frame_id_(frame_id),
    q_3d_ref_(Eigen::Vector3d::Zero()),
    q_3d_weight_(Eigen::Vector3d::Zero()),
    qf_3d_weight_(Eigen::Vector3d::Zero()),
    qi_3d_weight_(Eigen::Vector3d::Zero()) {
}


TaskSpace3DCost::TaskSpace3DCost()
  : CostFunctionComponentBase(),
    frame_id_(),
    q_3d_ref_(),
    q_3d_weight_(),
    qf_3d_weight_(),
    qi_3d_weight_() {
}


TaskSpace3DCost::~TaskSpace3DCost() {
}


void TaskSpace3DCost::set_q_3d_ref(const Eigen::Vector3d& q_3d_ref) {
  q_3d_ref_ = q_3d_ref;
}


void TaskSpace3DCost::set_q_weight(const Eigen::Vector3d& q_3d_weight) {
  q_3d_weight_ = q_3d_weight;
}


void TaskSpace3DCost::set_qf_weight(const Eigen::Vector3d& qf_3d_weight) {
  qf_3d_weight_ = qf_3d_weight;
}


void TaskSpace3DCost::set_qi_weight(const Eigen::Vector3d& qi_3d_weight) {
  qi_3d_weight_ = qi_3d_weight;
}


bool TaskSpace3DCost::useKinematics() const {
  return true;
}


double TaskSpace3DCost::evalStageCost(Robot& robot, 
                                      const ContactStatus& contact_status, 
                                      CostFunctionData& data, 
                                      const double t, const double dt, 
                                      const SplitSolution& s) const {
  double l = 0;
  data.diff_3d = robot.framePosition(frame_id_) - q_3d_ref_;
  l += (q_3d_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
  return 0.5 * dt * l;
}


void TaskSpace3DCost::evalStageCostDerivatives(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const double t, const double dt, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) const {
  data.J_6d.setZero();
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.J_3d.noalias() 
      = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
  kkt_residual.lq().noalias() 
      += dt * data.J_3d.transpose() * q_3d_weight_.asDiagonal() * data.diff_3d;
}


void TaskSpace3DCost::evalStageCostHessian(Robot& robot, 
                                           const ContactStatus& contact_status, 
                                           CostFunctionData& data, 
                                           const double t, const double dt, 
                                           const SplitSolution& s, 
                                           SplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qqq().noalias()
      += dt * data.J_3d.transpose() * q_3d_weight_.asDiagonal() * data.J_3d;
}


double TaskSpace3DCost::evalTerminalCost(Robot& robot, CostFunctionData& data, 
                                         const double t, 
                                         const SplitSolution& s) const {
  double l = 0;
  data.diff_3d = robot.framePosition(frame_id_) - q_3d_ref_;
  l += (qf_3d_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
  return 0.5 * l;
}


void TaskSpace3DCost::evalTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  data.J_6d.setZero();
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.J_3d.noalias() 
      = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
  kkt_residual.lq().noalias() 
      += data.J_3d.transpose() * qf_3d_weight_.asDiagonal() * data.diff_3d;
}


void TaskSpace3DCost::evalTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qqq().noalias()
      += data.J_3d.transpose() * qf_3d_weight_.asDiagonal() * data.J_3d;
}


double TaskSpace3DCost::evalImpulseCost(Robot& robot,  
                                        const ImpulseStatus& impulse_status,
                                        CostFunctionData& data, 
                                        const double t, 
                                        const ImpulseSplitSolution& s) const {
  double l = 0;
  data.diff_3d = robot.framePosition(frame_id_) - q_3d_ref_;
  l += (qi_3d_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
  return 0.5 * l;
}


void TaskSpace3DCost::evalImpulseCostDerivatives(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const double t, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  data.J_6d.setZero();
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.J_3d.noalias() 
      = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
  kkt_residual.lq().noalias() 
      += data.J_3d.transpose() * qi_3d_weight_.asDiagonal() * data.diff_3d;
}


void TaskSpace3DCost::evalImpulseCostHessian(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const double t, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qqq().noalias()
      += data.J_3d.transpose() * qi_3d_weight_.asDiagonal() * data.J_3d;
}

} // namespace robotoc