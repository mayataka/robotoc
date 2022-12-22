#include "robotoc/cost/task_space_3d_cost.hpp"


namespace robotoc {

TaskSpace3DCost::TaskSpace3DCost(const Robot& robot, const int frame_id)
  : CostFunctionComponentBase(),
    frame_id_(frame_id),
    const_ref_(Eigen::Vector3d::Zero()),
    weight_(Eigen::Vector3d::Zero()),
    weight_terminal_(Eigen::Vector3d::Zero()),
    weight_impact_(Eigen::Vector3d::Zero()),
    ref_(),
    use_nonconst_ref_(false),
    enable_cost_(false), 
    enable_cost_terminal_(false), 
    enable_cost_impact_(false) {
}


TaskSpace3DCost::TaskSpace3DCost(const Robot& robot, 
                                 const std::string& frame_name)
  : TaskSpace3DCost(robot, robot.frameId(frame_name)) {
}


TaskSpace3DCost::TaskSpace3DCost(const Robot& robot, const int frame_id,
                                 const std::shared_ptr<TaskSpace3DRefBase>& ref)
  : TaskSpace3DCost(robot, frame_id) {
  set_ref(ref);
}


TaskSpace3DCost::TaskSpace3DCost(const Robot& robot, const int frame_id,
                                 const Eigen::Vector3d& const_ref)
  : TaskSpace3DCost(robot, frame_id) {
  set_const_ref(const_ref);
}


TaskSpace3DCost::TaskSpace3DCost(const Robot& robot, 
                                 const std::string& frame_name,
                                 const std::shared_ptr<TaskSpace3DRefBase>& ref)
  : TaskSpace3DCost(robot, frame_name) {
  set_ref(ref);
}


TaskSpace3DCost::TaskSpace3DCost(const Robot& robot, 
                                 const std::string& frame_name,
                                 const Eigen::Vector3d& const_ref)
  : TaskSpace3DCost(robot, frame_name) {
  set_const_ref(const_ref);
}


TaskSpace3DCost::TaskSpace3DCost()
  : CostFunctionComponentBase(),
    frame_id_(0),
    const_ref_(Eigen::Vector3d::Zero()),
    weight_(Eigen::Vector3d::Zero()),
    weight_terminal_(Eigen::Vector3d::Zero()),
    weight_impact_(Eigen::Vector3d::Zero()),
    ref_(),
    use_nonconst_ref_(false),
    enable_cost_(false), 
    enable_cost_terminal_(false), 
    enable_cost_impact_(false) {
}


TaskSpace3DCost::~TaskSpace3DCost() {
}


void TaskSpace3DCost::set_ref(const std::shared_ptr<TaskSpace3DRefBase>& ref) {
  ref_ = ref;
  use_nonconst_ref_ = true;
}


void TaskSpace3DCost::set_const_ref(const Eigen::Vector3d& const_ref) {
  const_ref_ = const_ref;
  use_nonconst_ref_ = false;
}


void TaskSpace3DCost::set_weight(const Eigen::Vector3d& weight) {
  if (weight.minCoeff() < 0.0) {
    throw std::invalid_argument(
        "[TaskSpace3DCost] invalid argument: elements of 'weight' must be non-negative!");
  }
  weight_ = weight;
  enable_cost_ = (!weight.isZero());
}


void TaskSpace3DCost::set_weight_terminal(const Eigen::Vector3d& weight_terminal) {
  if (weight_terminal.minCoeff() < 0.0) {
    throw std::invalid_argument(
        "[TaskSpace3DCost] invalid argument: elements of 'weight_terminal' must be non-negative!");
  }
  weight_terminal_ = weight_terminal;
  enable_cost_terminal_ = (!weight_terminal.isZero());
}


void TaskSpace3DCost::set_weight_impact(const Eigen::Vector3d& weight_impact) {
  if (weight_impact.minCoeff() < 0.0) {
    throw std::invalid_argument(
        "[TaskSpace3DCost] invalid argument: elements of 'weight_impact' must be non-negative!");
  }
  weight_impact_ = weight_impact;
  enable_cost_impact_ = (!weight_impact.isZero());
}


double TaskSpace3DCost::evalStageCost(Robot& robot, 
                                      const ContactStatus& contact_status, 
                                      const GridInfo& grid_info, 
                                      const SplitSolution& s,
                                      CostFunctionData& data) const {
  if (enable_cost_ && isCostActive(grid_info)) {
    evalDiff(robot, data, grid_info);
    const double l = (weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
    return 0.5 * grid_info.dt * l;
  }
  else {
    return 0.0;
  }
}


void TaskSpace3DCost::evalStageCostDerivatives(
    Robot& robot, const ContactStatus& contact_status, const GridInfo& grid_info,
    const SplitSolution& s, CostFunctionData& data, 
    SplitKKTResidual& kkt_residual) const {
  if (enable_cost_ && isCostActive(grid_info)) {
    data.J_6d.setZero();
    robot.getFrameJacobian(frame_id_, data.J_6d);
    data.J_3d.noalias() 
        = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
    kkt_residual.lq().noalias() 
        += grid_info.dt * data.J_3d.transpose() * weight_.asDiagonal() * data.diff_3d;
  }
}


void TaskSpace3DCost::evalStageCostHessian(Robot& robot, 
                                           const ContactStatus& contact_status, 
                                           const GridInfo& grid_info, 
                                           const SplitSolution& s, 
                                           CostFunctionData& data, 
                                           SplitKKTMatrix& kkt_matrix) const {
  if (enable_cost_ && isCostActive(grid_info)) {
    kkt_matrix.Qqq().noalias()
        += grid_info.dt * data.J_3d.transpose() * weight_.asDiagonal() * data.J_3d;
  }
}


double TaskSpace3DCost::evalTerminalCost(Robot& robot, const GridInfo& grid_info, 
                                         const SplitSolution& s,
                                         CostFunctionData& data) const {
  if (enable_cost_terminal_ && isCostActive(grid_info)) {
    evalDiff(robot, data, grid_info);
    const double l = (weight_terminal_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
    return 0.5 * l;
  }
  else {
    return 0.0;
  }
}


void TaskSpace3DCost::evalTerminalCostDerivatives(Robot& robot, const GridInfo& grid_info, 
                                                  const SplitSolution& s, CostFunctionData& data,
                                                  SplitKKTResidual& kkt_residual) const {
  if (enable_cost_terminal_ && isCostActive(grid_info)) {
    data.J_6d.setZero();
    robot.getFrameJacobian(frame_id_, data.J_6d);
    data.J_3d.noalias() 
        = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
    kkt_residual.lq().noalias() 
        += data.J_3d.transpose() * weight_terminal_.asDiagonal() * data.diff_3d;
  }
}


void TaskSpace3DCost::evalTerminalCostHessian(Robot& robot, const GridInfo& grid_info,
                                              const SplitSolution& s, CostFunctionData& data, 
                                              SplitKKTMatrix& kkt_matrix) const {
  if (enable_cost_terminal_ && isCostActive(grid_info)) {
    kkt_matrix.Qqq().noalias()
        += data.J_3d.transpose() * weight_terminal_.asDiagonal() * data.J_3d;
  }
}


double TaskSpace3DCost::evalImpactCost(Robot& robot, const ImpactStatus& impact_status,
                                       const GridInfo& grid_info, const SplitSolution& s,
                                       CostFunctionData& data) const {
  if (enable_cost_impact_ && isCostActive(grid_info)) {
    evalDiff(robot, data, grid_info);
    const double l = (weight_impact_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
    return 0.5 * l;
  }
  else {
    return 0.0;
  }
}


void TaskSpace3DCost::evalImpactCostDerivatives(
    Robot& robot, const ImpactStatus& impact_status, const GridInfo& grid_info,
     const SplitSolution& s, CostFunctionData& data, 
    SplitKKTResidual& kkt_residual) const {
  if (enable_cost_impact_ && isCostActive(grid_info)) {
    data.J_6d.setZero();
    robot.getFrameJacobian(frame_id_, data.J_6d);
    data.J_3d.noalias() 
        = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
    kkt_residual.lq().noalias() 
        += data.J_3d.transpose() * weight_impact_.asDiagonal() * data.diff_3d;
  }
}


void TaskSpace3DCost::evalImpactCostHessian(
    Robot& robot, const ImpactStatus& impact_status, const GridInfo& grid_info,
    const SplitSolution& s, CostFunctionData& data, 
    SplitKKTMatrix& kkt_matrix) const {
  if (enable_cost_impact_ && isCostActive(grid_info)) {
    kkt_matrix.Qqq().noalias()
        += data.J_3d.transpose() * weight_impact_.asDiagonal() * data.J_3d;
  }
}

} // namespace robotoc