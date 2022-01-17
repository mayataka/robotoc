#include "robotoc/cost/task_space_6d_cost.hpp"


namespace robotoc {

TaskSpace6DCost::TaskSpace6DCost(const Robot& robot, const int frame_id)
  : CostFunctionComponentBase(),
    frame_id_(frame_id),
    x6d_ref_(SE3(Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero())),
    x6d_ref_inv_(x6d_ref_.inverse()),
    x6d_weight_(Eigen::VectorXd::Zero(6)), 
    x6df_weight_(Eigen::VectorXd::Zero(6)), 
    x6di_weight_(Eigen::VectorXd::Zero(6)) {
}


TaskSpace6DCost::TaskSpace6DCost()
  : CostFunctionComponentBase(),
    frame_id_(0),
    x6d_ref_(SE3(Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero())),
    x6d_ref_inv_(x6d_ref_.inverse()),
    x6d_weight_(Eigen::VectorXd::Zero(6)), 
    x6df_weight_(Eigen::VectorXd::Zero(6)), 
    x6di_weight_(Eigen::VectorXd::Zero(6)) {
}


TaskSpace6DCost::~TaskSpace6DCost() {
}


void TaskSpace6DCost::set_x6d_ref(const Eigen::Vector3d& trans_ref, 
                                  const Eigen::Matrix3d& rot_ref) {
  x6d_ref_ = SE3(rot_ref, trans_ref);
  x6d_ref_inv_ = x6d_ref_.inverse();
}


void TaskSpace6DCost::set_x6d_weight(const Eigen::Vector3d& trans_weight, 
                                     const Eigen::Vector3d& rot_weight) {
  x6d_weight_.template head<3>() = rot_weight;
  x6d_weight_.template tail<3>() = trans_weight;
}


void TaskSpace6DCost::set_x6df_weight(const Eigen::Vector3d& trans_weight, 
                                      const Eigen::Vector3d& rot_weight) {
  x6df_weight_.template head<3>() = rot_weight;
  x6df_weight_.template tail<3>() = trans_weight;
}


void TaskSpace6DCost::set_x6di_weight(const Eigen::Vector3d& trans_weight, 
                                      const Eigen::Vector3d& rot_weight) {
  x6di_weight_.template head<3>() = rot_weight;
  x6di_weight_.template tail<3>() = trans_weight;
}
 

bool TaskSpace6DCost::useKinematics() const {
  return true;
}


double TaskSpace6DCost::evalStageCost(Robot& robot, 
                                      const ContactStatus& contact_status, 
                                      CostFunctionData& data, 
                                      const GridInfo& grid_info, 
                                      const SplitSolution& s) const {
  double l = 0;
  data.diff_x6d = x6d_ref_inv_ * robot.framePlacement(frame_id_);
  data.diff_6d = Log6Map(data.diff_x6d);
  l += (x6d_weight_.array()*data.diff_6d.array()*data.diff_6d.array()).sum();
  return 0.5 * grid_info.dt * l;
}


void TaskSpace6DCost::evalStageCostDerivatives(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const GridInfo& grid_info, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) const {
  data.J_66.setZero();
  computeJLog6Map(data.diff_x6d, data.J_66);
  data.J_6d.setZero();
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.JJ_6d.noalias() = data.J_66 * data.J_6d;
  kkt_residual.lq().noalias() 
      += grid_info.dt * data.JJ_6d.transpose() * x6d_weight_.asDiagonal() * data.diff_6d;
}


void TaskSpace6DCost::evalStageCostHessian(Robot& robot, 
                                           const ContactStatus& contact_status, 
                                           CostFunctionData& data, 
                                           const GridInfo& grid_info, 
                                           const SplitSolution& s, 
                                           SplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qqq().noalias()
      += grid_info.dt * data.JJ_6d.transpose() * x6d_weight_.asDiagonal() * data.JJ_6d;
}


double TaskSpace6DCost::evalTerminalCost(Robot& robot, CostFunctionData& data, 
                                         const GridInfo& grid_info, 
                                         const SplitSolution& s) const {
  double l = 0;
  data.diff_x6d = x6d_ref_inv_ * robot.framePlacement(frame_id_);
  data.diff_6d = Log6Map(data.diff_x6d);
  l += (x6df_weight_.array()*data.diff_6d.array()*data.diff_6d.array()).sum();
  return 0.5 * l;
}


void TaskSpace6DCost::evalTerminalCostDerivatives(
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


void TaskSpace6DCost::evalTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const GridInfo& grid_info, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qqq().noalias()
      += data.JJ_6d.transpose() * x6df_weight_.asDiagonal() * data.JJ_6d;
}


double TaskSpace6DCost::evalImpulseCost(Robot& robot, 
                                        const ImpulseStatus& impulse_status,
                                        CostFunctionData& data, 
                                        const GridInfo& grid_info, 
                                        const ImpulseSplitSolution& s) const {
  double l = 0;
  data.diff_x6d = x6d_ref_inv_ * robot.framePlacement(frame_id_);
  data.diff_6d = Log6Map(data.diff_x6d);
  l += (x6di_weight_.array()*data.diff_6d.array()*data.diff_6d.array()).sum();
  return 0.5 * l;
}


void TaskSpace6DCost::evalImpulseCostDerivatives(
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


void TaskSpace6DCost::evalImpulseCostHessian(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qqq().noalias()
      += data.JJ_6d.transpose() * x6di_weight_.asDiagonal() * data.JJ_6d;
}

} // namespace robotoc