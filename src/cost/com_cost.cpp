#include "robotoc/cost/com_cost.hpp"


namespace robotoc {

CoMCost::CoMCost(const Robot& robot)
  : CostFunctionComponentBase(),
    const_ref_(Eigen::Vector3d::Zero()),
    weight_(Eigen::Vector3d::Zero()),
    weight_terminal_(Eigen::Vector3d::Zero()),
    weight_impulse_(Eigen::Vector3d::Zero()),
    ref_(),
    use_nonconst_ref_(false),
    enable_cost_(false), 
    enable_cost_terminal_(false), 
    enable_cost_impulse_(false) {
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
    weight_impulse_(Eigen::Vector3d::Zero()),
    ref_(),
    use_nonconst_ref_(false),
    enable_cost_(false), 
    enable_cost_terminal_(false), 
    enable_cost_impulse_(false) {
}


CoMCost::~CoMCost() {
}


std::shared_ptr<CostFunctionComponentBase> CoMCost::clone() const {
  auto cost = std::make_shared<CoMCost>(*this);
  if (use_nonconst_ref_) {
    auto ref = ref_->clone();
    cost->set_ref(ref);
  }
  return cost;
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
  try {
    if (weight.minCoeff() < 0.0) {
      throw std::invalid_argument(
          "invalid argument: elements of weight must be non-negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  weight_ = weight;
  enable_cost_ = (!weight.isZero());
}


void CoMCost::set_weight_terminal(const Eigen::Vector3d& weight_terminal) {
  try {
    if (weight_terminal.minCoeff() < 0.0) {
      throw std::invalid_argument(
          "invalid argument: elements of weight must be non-negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  weight_terminal_ = weight_terminal;
  enable_cost_terminal_ = (!weight_terminal.isZero());
}


void CoMCost::set_weight_impulse(const Eigen::Vector3d& weight_impulse) {
  try {
    if (weight_impulse.minCoeff() < 0.0) {
      throw std::invalid_argument(
          "invalid argument: elements of weight must be non-negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  weight_impulse_ = weight_impulse;
  enable_cost_impulse_ = (!weight_impulse.isZero());
}


bool CoMCost::useKinematics() const {
  return true;
}


double CoMCost::evalStageCost(Robot& robot, const ContactStatus& contact_status, 
                              CostFunctionData& data, const GridInfo& grid_info, 
                              const SplitSolution& s) const {
  if (enable_cost_ && isCostActive(grid_info)) {
    evalDiff(robot, data, grid_info);
    const double l = (weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
    return 0.5 * grid_info.dt * l;
  }
  else {
    return 0.0;
  }
}


void CoMCost::evalStageCostDerivatives(Robot& robot, 
                                       const ContactStatus& contact_status, 
                                       CostFunctionData& data, 
                                       const GridInfo& grid_info,
                                       const SplitSolution& s,
                                       SplitKKTResidual& kkt_residual) const {
  if (enable_cost_ && isCostActive(grid_info)) {
    data.J_3d.setZero();
    robot.getCoMJacobian(data.J_3d);
    kkt_residual.lq().noalias() 
        += grid_info.dt * data.J_3d.transpose() * weight_.asDiagonal() * data.diff_3d;
  }
}


void CoMCost::evalStageCostHessian(Robot& robot, 
                                   const ContactStatus& contact_status, 
                                   CostFunctionData& data, 
                                   const GridInfo& grid_info,
                                   const SplitSolution& s, 
                                   SplitKKTMatrix& kkt_matrix) const {
  if (enable_cost_ && isCostActive(grid_info)) {
    kkt_matrix.Qqq().noalias()
        += grid_info.dt * data.J_3d.transpose() * weight_.asDiagonal() * data.J_3d;
  }
}


double CoMCost::evalTerminalCost(Robot& robot, CostFunctionData& data, 
                                 const GridInfo& grid_info,
                                 const SplitSolution& s) const {
  if (enable_cost_terminal_ && isCostActive(grid_info)) {
    evalDiff(robot, data, grid_info);
    const double l = (weight_terminal_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
    return 0.5 * l;
  }
  else {
    return 0.0;
  }
}


void CoMCost::evalTerminalCostDerivatives(Robot& robot, CostFunctionData& data, 
                                          const GridInfo& grid_info,
                                          const SplitSolution& s, 
                                          SplitKKTResidual& kkt_residual) const {
  if (enable_cost_terminal_ && isCostActive(grid_info)) {
    data.J_3d.setZero();
    robot.getCoMJacobian(data.J_3d);
    kkt_residual.lq().noalias() 
        += data.J_3d.transpose() * weight_terminal_.asDiagonal() * data.diff_3d;
  }
}


void CoMCost::evalTerminalCostHessian(Robot& robot, CostFunctionData& data, 
                                      const GridInfo& grid_info,
                                      const SplitSolution& s, 
                                      SplitKKTMatrix& kkt_matrix) const {
  if (enable_cost_terminal_ && isCostActive(grid_info)) {
    data.J_3d.setZero();
    robot.getCoMJacobian(data.J_3d);
    kkt_matrix.Qqq().noalias()
        += data.J_3d.transpose() * weight_terminal_.asDiagonal() * data.J_3d;
  }
}


double CoMCost::evalImpulseCost(Robot& robot, const ImpulseStatus& impulse_status, 
                                CostFunctionData& data, 
                                const GridInfo& grid_info,
                                const ImpulseSplitSolution& s) const {
  if (enable_cost_impulse_ && isCostActive(grid_info)) {
    evalDiff(robot, data, grid_info);
    const double l = (weight_impulse_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
    return 0.5 * l;
  }
  else {
    return 0.0;
  }
}


void CoMCost::evalImpulseCostDerivatives(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  if (enable_cost_impulse_ && isCostActive(grid_info)) {
    data.J_3d.setZero();
    robot.getCoMJacobian(data.J_3d);
    kkt_residual.lq().noalias() 
        += data.J_3d.transpose() * weight_impulse_.asDiagonal() * data.diff_3d;
  }
}


void CoMCost::evalImpulseCostHessian(Robot& robot, 
                                     const ImpulseStatus& impulse_status, 
                                     CostFunctionData& data, 
                                     const GridInfo& grid_info,
                                     const ImpulseSplitSolution& s, 
                                     ImpulseSplitKKTMatrix& kkt_matrix) const {
  if (enable_cost_impulse_ && isCostActive(grid_info)) {
    kkt_matrix.Qqq().noalias()
        += data.J_3d.transpose() * weight_impulse_.asDiagonal() * data.J_3d;
  }
}

} // namespace robotoc