#include "robotoc/cost/task_space_3d_cost.hpp"


namespace robotoc {

TaskSpace3DCost::TaskSpace3DCost(const Robot& robot, const int frame_id)
  : CostFunctionComponentBase(),
    frame_id_(frame_id),
    x3d_ref_(Eigen::Vector3d::Zero()),
    x3d_weight_(Eigen::Vector3d::Zero()),
    x3df_weight_(Eigen::Vector3d::Zero()),
    x3di_weight_(Eigen::Vector3d::Zero()) {
}


TaskSpace3DCost::TaskSpace3DCost()
  : CostFunctionComponentBase(),
    frame_id_(),
    x3d_ref_(),
    x3d_weight_(),
    x3df_weight_(),
    x3di_weight_() {
}


TaskSpace3DCost::~TaskSpace3DCost() {
}


void TaskSpace3DCost::set_x3d_ref(const Eigen::Vector3d& x3d_ref) {
  x3d_ref_ = x3d_ref;
}


void TaskSpace3DCost::set_x3d_weight(const Eigen::Vector3d& x3d_weight) {
  x3d_weight_ = x3d_weight;
}


void TaskSpace3DCost::set_x3df_weight(const Eigen::Vector3d& x3df_weight) {
  x3df_weight_ = x3df_weight;
}


void TaskSpace3DCost::set_x3di_weight(const Eigen::Vector3d& x3di_weight) {
  x3di_weight_ = x3di_weight;
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
  data.diff_3d = robot.framePosition(frame_id_) - x3d_ref_;
  l += (x3d_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
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
      += dt * data.J_3d.transpose() * x3d_weight_.asDiagonal() * data.diff_3d;
}


void TaskSpace3DCost::evalStageCostHessian(Robot& robot, 
                                           const ContactStatus& contact_status, 
                                           CostFunctionData& data, 
                                           const double t, const double dt, 
                                           const SplitSolution& s, 
                                           SplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qqq().noalias()
      += dt * data.J_3d.transpose() * x3d_weight_.asDiagonal() * data.J_3d;
}


double TaskSpace3DCost::evalTerminalCost(Robot& robot, CostFunctionData& data, 
                                         const double t, 
                                         const SplitSolution& s) const {
  double l = 0;
  data.diff_3d = robot.framePosition(frame_id_) - x3d_ref_;
  l += (x3df_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
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
      += data.J_3d.transpose() * x3df_weight_.asDiagonal() * data.diff_3d;
}


void TaskSpace3DCost::evalTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qqq().noalias()
      += data.J_3d.transpose() * x3df_weight_.asDiagonal() * data.J_3d;
}


double TaskSpace3DCost::evalImpulseCost(Robot& robot,  
                                        const ImpulseStatus& impulse_status,
                                        CostFunctionData& data, 
                                        const double t, 
                                        const ImpulseSplitSolution& s) const {
  double l = 0;
  data.diff_3d = robot.framePosition(frame_id_) - x3d_ref_;
  l += (x3di_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
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
      += data.J_3d.transpose() * x3di_weight_.asDiagonal() * data.diff_3d;
}


void TaskSpace3DCost::evalImpulseCostHessian(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const double t, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qqq().noalias()
      += data.J_3d.transpose() * x3di_weight_.asDiagonal() * data.J_3d;
}

} // namespace robotoc