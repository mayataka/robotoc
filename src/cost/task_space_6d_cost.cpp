#include "robotoc/cost/task_space_6d_cost.hpp"


namespace robotoc {

TaskSpace6DCost::TaskSpace6DCost(const Robot& robot, const int frame_id)
  : CostFunctionComponentBase(),
    frame_id_(frame_id),
    const_ref_(SE3(Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero())),
    const_ref_inv_(const_ref_.inverse()),
    weight_(Eigen::VectorXd::Zero(6)), 
    weight_terminal_(Eigen::VectorXd::Zero(6)), 
    weight_impulse_(Eigen::VectorXd::Zero(6)),
    ref_(),
    use_nonconst_ref_(false) {
}


TaskSpace6DCost::TaskSpace6DCost(const Robot& robot, 
                                 const std::string& frame_name)
  : TaskSpace6DCost(robot, robot.frameId(frame_name)) {
}


TaskSpace6DCost::TaskSpace6DCost(const Robot& robot, const int frame_id,
                                 const std::shared_ptr<TaskSpace6DRefBase>& ref)
  : TaskSpace6DCost(robot, frame_id) {
  set_ref(ref);
}


TaskSpace6DCost::TaskSpace6DCost(const Robot& robot, const int frame_id,
                                 const SE3& const_ref)
  : TaskSpace6DCost(robot, frame_id) {
  set_const_ref(const_ref);
}


TaskSpace6DCost::TaskSpace6DCost(const Robot& robot, const int frame_id,
                                 const Eigen::Vector3d& const_position_ref,
                                 const Eigen::Matrix3d& const_rotation_ref) 
  : TaskSpace6DCost(robot, frame_id) {
  set_const_ref(const_position_ref, const_rotation_ref);
}


TaskSpace6DCost::TaskSpace6DCost(const Robot& robot, 
                                 const std::string& frame_name,
                                 const std::shared_ptr<TaskSpace6DRefBase>& ref)
  : TaskSpace6DCost(robot, frame_name) {
  set_ref(ref);
}


TaskSpace6DCost::TaskSpace6DCost(const Robot& robot,
                                 const std::string& frame_name,
                                 const SE3& const_ref)
  : TaskSpace6DCost(robot, frame_name) {
  set_const_ref(const_ref);
}


TaskSpace6DCost::TaskSpace6DCost(const Robot& robot,
                                 const std::string& frame_name,
                                 const Eigen::Vector3d& const_position_ref,
                                 const Eigen::Matrix3d& const_rotation_ref) 
  : TaskSpace6DCost(robot, frame_name) {
  set_const_ref(const_position_ref, const_rotation_ref);
}


TaskSpace6DCost::TaskSpace6DCost()
  : CostFunctionComponentBase(),
    frame_id_(0),
    const_ref_(SE3(Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero())),
    const_ref_inv_(const_ref_.inverse()),
    weight_(Eigen::VectorXd::Zero(6)), 
    weight_terminal_(Eigen::VectorXd::Zero(6)), 
    weight_impulse_(Eigen::VectorXd::Zero(6)),
    ref_(),
    use_nonconst_ref_(false) {
}


TaskSpace6DCost::~TaskSpace6DCost() {
}


void TaskSpace6DCost::set_ref(const std::shared_ptr<TaskSpace6DRefBase>& ref) {
  ref_ = ref;
  use_nonconst_ref_ = true;
}


void TaskSpace6DCost::set_const_ref(const SE3& const_ref) {
  const_ref_ = const_ref;
  const_ref_inv_ = const_ref.inverse();
  use_nonconst_ref_ = false;
}


void TaskSpace6DCost::set_const_ref(const Eigen::Vector3d& const_position_ref,
                                    const Eigen::Matrix3d& const_rotation_ref) {
  set_const_ref(SE3(const_rotation_ref, const_position_ref));
}


void TaskSpace6DCost::set_weight(const Eigen::Vector3d& weight_position,
                                 const Eigen::Vector3d& weight_rotation) {
  weight_.template head<3>() = weight_rotation;
  weight_.template tail<3>() = weight_position;
}


void TaskSpace6DCost::set_weight_terminal(
    const Eigen::Vector3d& weight_position_terminal, 
    const Eigen::Vector3d& weight_rotation_terminal) {
  weight_terminal_.template head<3>() = weight_rotation_terminal;
  weight_terminal_.template tail<3>() = weight_position_terminal;
}


void TaskSpace6DCost::set_weight_impulse(
    const Eigen::Vector3d& weight_position_impulse, 
    const Eigen::Vector3d& weight_rotation_impulse) {
  weight_impulse_.template head<3>() = weight_rotation_impulse;
  weight_impulse_.template tail<3>() = weight_position_impulse;
}


bool TaskSpace6DCost::useKinematics() const {
  return true;
}


double TaskSpace6DCost::evalStageCost(Robot& robot, 
                                      const ContactStatus& contact_status, 
                                      CostFunctionData& data, 
                                      const GridInfo& grid_info, 
                                      const SplitSolution& s) const {
  if (isCostActive(grid_info)) {
    evalDiff(robot, data, grid_info);
    const double l = (weight_.array()*data.diff_6d.array()*data.diff_6d.array()).sum();
    return 0.5 * grid_info.dt * l;
  }
  else {
    return 0.0;
  }
}


void TaskSpace6DCost::evalStageCostDerivatives(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const GridInfo& grid_info, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) const {
  if (isCostActive(grid_info)) {
    data.J_66.setZero();
    computeJLog6Map(data.diff_x6d, data.J_66);
    data.J_6d.setZero();
    robot.getFrameJacobian(frame_id_, data.J_6d);
    data.JJ_6d.noalias() = data.J_66 * data.J_6d;
    kkt_residual.lq().noalias() 
        += grid_info.dt * data.JJ_6d.transpose() * weight_.asDiagonal() * data.diff_6d;
  }
}


void TaskSpace6DCost::evalStageCostHessian(Robot& robot, 
                                           const ContactStatus& contact_status, 
                                           CostFunctionData& data, 
                                           const GridInfo& grid_info, 
                                           const SplitSolution& s, 
                                           SplitKKTMatrix& kkt_matrix) const {
  if (isCostActive(grid_info)) {
    kkt_matrix.Qqq().noalias()
        += grid_info.dt * data.JJ_6d.transpose() * weight_.asDiagonal() * data.JJ_6d;
  }
}


double TaskSpace6DCost::evalTerminalCost(Robot& robot, CostFunctionData& data, 
                                         const GridInfo& grid_info, 
                                         const SplitSolution& s) const {
  if (isCostActive(grid_info)) {
    evalDiff(robot, data, grid_info);
    const double l = (weight_terminal_.array()*data.diff_6d.array()*data.diff_6d.array()).sum();
    return 0.5 * l;
  }
  else {
    return 0.0;
  }
}


void TaskSpace6DCost::evalTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const GridInfo& grid_info, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  if (isCostActive(grid_info)) {
    data.J_66.setZero();
    computeJLog6Map(data.diff_x6d, data.J_66);
    data.J_6d.setZero();
    robot.getFrameJacobian(frame_id_, data.J_6d);
    data.JJ_6d.noalias() = data.J_66 * data.J_6d;
    kkt_residual.lq().noalias() 
        += data.JJ_6d.transpose() * weight_terminal_.asDiagonal() * data.diff_6d;
  }
}


void TaskSpace6DCost::evalTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const GridInfo& grid_info, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  if (isCostActive(grid_info)) {
    kkt_matrix.Qqq().noalias()
        += data.JJ_6d.transpose() * weight_terminal_.asDiagonal() * data.JJ_6d;
  }
}


double TaskSpace6DCost::evalImpulseCost(Robot& robot, 
                                        const ImpulseStatus& impulse_status,
                                        CostFunctionData& data, 
                                        const GridInfo& grid_info, 
                                        const ImpulseSplitSolution& s) const {
  if (isCostActive(grid_info)) {
    evalDiff(robot, data, grid_info);
    const double l = (weight_impulse_.array()*data.diff_6d.array()*data.diff_6d.array()).sum();
    return 0.5 * l;
  }
  else {
    return 0.0;
  }
}


void TaskSpace6DCost::evalImpulseCostDerivatives(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  if (isCostActive(grid_info)) {
    data.J_66.setZero();
    computeJLog6Map(data.diff_x6d, data.J_66);
    data.J_6d.setZero();
    robot.getFrameJacobian(frame_id_, data.J_6d);
    data.JJ_6d.noalias() = data.J_66 * data.J_6d;
    kkt_residual.lq().noalias() 
        += data.JJ_6d.transpose() * weight_impulse_.asDiagonal() * data.diff_6d;
  }
}


void TaskSpace6DCost::evalImpulseCostHessian(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTMatrix& kkt_matrix) const {
  if (isCostActive(grid_info)) {
    kkt_matrix.Qqq().noalias()
        += data.JJ_6d.transpose() * weight_impulse_.asDiagonal() * data.JJ_6d;
  }
}

} // namespace robotoc