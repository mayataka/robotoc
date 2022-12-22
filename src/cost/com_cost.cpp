#include "robotoc/cost/com_cost.hpp"


namespace robotoc {

CoMCost::CoMCost(const Robot& robot)
  : CostFunctionComponentBase(),
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


CoMCost::CoMCost(const Robot& robot, const std::shared_ptr<CoMRefBase>& ref)
  : CoMCost(robot) {
  set_ref(ref);
}


CoMCost::CoMCost(const Robot& robot, const Eigen::Vector3d& const_ref)
  : CoMCost(robot) {
  set_const_ref(const_ref);
}


CoMCost::CoMCost()
  : CostFunctionComponentBase(),
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


CoMCost::~CoMCost() {
}


void CoMCost::set_ref(const std::shared_ptr<CoMRefBase>& ref) {
  ref_ = ref;
  use_nonconst_ref_ = true;
}


void CoMCost::set_const_ref(const Eigen::Vector3d& const_ref) {
  const_ref_ = const_ref;
  use_nonconst_ref_ = false;
}


void CoMCost::set_weight(const Eigen::Vector3d& weight) {
  if (weight.minCoeff() < 0.0) {
    throw std::invalid_argument(
        "[CoMCost] invalid argument: elements of 'weight' must be non-negative!");
  }
  weight_ = weight;
  enable_cost_ = (!weight.isZero());
}


void CoMCost::set_weight_terminal(const Eigen::Vector3d& weight_terminal) {
  if (weight_terminal.minCoeff() < 0.0) {
    throw std::invalid_argument(
        "[CoMCost] invalid argument: elements of 'weight_terminal' must be non-negative!");
  }
  weight_terminal_ = weight_terminal;
  enable_cost_terminal_ = (!weight_terminal.isZero());
}


void CoMCost::set_weight_impact(const Eigen::Vector3d& weight_impact) {
  if (weight_impact.minCoeff() < 0.0) {
    throw std::invalid_argument(
        "[CoMCost] invalid argument: elements of 'weight_impact' must be non-negative!");
  }
  weight_impact_ = weight_impact;
  enable_cost_impact_ = (!weight_impact.isZero());
}


double CoMCost::evalStageCost(Robot& robot, const ContactStatus& contact_status, 
                              const GridInfo& grid_info, const SplitSolution& s,
                              CostFunctionData& data) const {
  if (enable_cost_ && isCostActive(grid_info)) {
    evalDiff(robot, data, grid_info);
    const double l = (weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
    return 0.5 * l;
  }
  else {
    return 0.0;
  }
}


void CoMCost::evalStageCostDerivatives(Robot& robot, 
                                       const ContactStatus& contact_status, 
                                       const GridInfo& grid_info,
                                       const SplitSolution& s,
                                       CostFunctionData& data, 
                                       SplitKKTResidual& kkt_residual) const {
  if (enable_cost_ && isCostActive(grid_info)) {
    data.J_3d.setZero();
    robot.getCoMJacobian(data.J_3d);
    kkt_residual.lq().noalias() 
        += data.J_3d.transpose() * weight_.asDiagonal() * data.diff_3d;
  }
}


void CoMCost::evalStageCostHessian(Robot& robot, 
                                   const ContactStatus& contact_status, 
                                   const GridInfo& grid_info,
                                   const SplitSolution& s, 
                                   CostFunctionData& data, 
                                   SplitKKTMatrix& kkt_matrix) const {
  if (enable_cost_ && isCostActive(grid_info)) {
    kkt_matrix.Qqq().noalias()
        += data.J_3d.transpose() * weight_.asDiagonal() * data.J_3d;
  }
}


double CoMCost::evalTerminalCost(Robot& robot, const GridInfo& grid_info,
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


void CoMCost::evalTerminalCostDerivatives(Robot& robot, const GridInfo& grid_info,
                                          const SplitSolution& s,
                                          CostFunctionData& data, 
                                          SplitKKTResidual& kkt_residual) const {
  if (enable_cost_terminal_ && isCostActive(grid_info)) {
    data.J_3d.setZero();
    robot.getCoMJacobian(data.J_3d);
    kkt_residual.lq().noalias() 
        += data.J_3d.transpose() * weight_terminal_.asDiagonal() * data.diff_3d;
  }
}


void CoMCost::evalTerminalCostHessian(Robot& robot, const GridInfo& grid_info,
                                      const SplitSolution& s, 
                                      CostFunctionData& data, 
                                      SplitKKTMatrix& kkt_matrix) const {
  if (enable_cost_terminal_ && isCostActive(grid_info)) {
    data.J_3d.setZero();
    robot.getCoMJacobian(data.J_3d);
    kkt_matrix.Qqq().noalias()
        += data.J_3d.transpose() * weight_terminal_.asDiagonal() * data.J_3d;
  }
}


double CoMCost::evalImpactCost(Robot& robot, const ImpactStatus& impact_status, 
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


void CoMCost::evalImpactCostDerivatives(Robot& robot,
                                        const ImpactStatus& impact_status,
                                        const GridInfo& grid_info,
                                        const SplitSolution& s,
                                        CostFunctionData& data,
                                        SplitKKTResidual& kkt_residual) const {
  if (enable_cost_impact_ && isCostActive(grid_info)) {
    data.J_3d.setZero();
    robot.getCoMJacobian(data.J_3d);
    kkt_residual.lq().noalias() 
        += data.J_3d.transpose() * weight_impact_.asDiagonal() * data.diff_3d;
  }
}


void CoMCost::evalImpactCostHessian(Robot& robot, 
                                    const ImpactStatus& impact_status, 
                                    const GridInfo& grid_info,
                                    const SplitSolution& s, 
                                    CostFunctionData& data, 
                                    SplitKKTMatrix& kkt_matrix) const {
  if (enable_cost_impact_ && isCostActive(grid_info)) {
    kkt_matrix.Qqq().noalias()
        += data.J_3d.transpose() * weight_impact_.asDiagonal() * data.J_3d;
  }
}

} // namespace robotoc