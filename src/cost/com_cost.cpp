#include "robotoc/cost/com_cost.hpp"


namespace robotoc {

CoMCost::CoMCost(const Robot& robot)
  : CostFunctionComponentBase(),
    com_ref_(Eigen::Vector3d::Zero()),
    com_weight_(Eigen::Vector3d::Zero()),
    comf_weight_(Eigen::Vector3d::Zero()),
    comi_weight_(Eigen::Vector3d::Zero()) {
}


CoMCost::CoMCost()
  : CostFunctionComponentBase(),
    com_ref_(),
    com_weight_(),
    comf_weight_(),
    comi_weight_() {
}


CoMCost::~CoMCost() {
}


void CoMCost::set_com_ref(const Eigen::Vector3d& com_ref) {
  com_ref_ = com_ref;
}


void CoMCost::set_com_weight(const Eigen::Vector3d& com_weight) {
  com_weight_ = com_weight;
}


void CoMCost::set_comf_weight(const Eigen::Vector3d& comf_weight) {
  comf_weight_ = comf_weight;
}


void CoMCost::set_comi_weight(const Eigen::Vector3d& comi_weight) {
  comi_weight_ = comi_weight;
}


bool CoMCost::useKinematics() const {
  return true;
}


double CoMCost::evalStageCost(Robot& robot, const ContactStatus& contact_status, 
                              CostFunctionData& data, const GridInfo& grid_info, 
                              const SplitSolution& s) const {
  double l = 0;
  data.diff_3d = robot.CoM() - com_ref_;
  l += (com_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
  return 0.5 * grid_info.dt * l;
}


void CoMCost::evalStageCostDerivatives(Robot& robot, 
                                       const ContactStatus& contact_status, 
                                       CostFunctionData& data, 
                                       const GridInfo& grid_info,
                                       const SplitSolution& s,
                                       SplitKKTResidual& kkt_residual) const {
  data.J_3d.setZero();
  robot.getCoMJacobian(data.J_3d);
  kkt_residual.lq().noalias() 
      += grid_info.dt * data.J_3d.transpose() * com_weight_.asDiagonal() * data.diff_3d;
}


void CoMCost::evalStageCostHessian(Robot& robot, 
                                   const ContactStatus& contact_status, 
                                   CostFunctionData& data, 
                                   const GridInfo& grid_info,
                                   const SplitSolution& s, 
                                   SplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qqq().noalias()
      += grid_info.dt * data.J_3d.transpose() * com_weight_.asDiagonal() * data.J_3d;
}


double CoMCost::evalTerminalCost(Robot& robot, CostFunctionData& data, 
                                 const GridInfo& grid_info,
                                 const SplitSolution& s) const {
  double l = 0;
  data.diff_3d = robot.CoM() - com_ref_;
  l += (comf_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
  return 0.5 * l;
}


void CoMCost::evalTerminalCostDerivatives(Robot& robot, CostFunctionData& data, 
                                          const GridInfo& grid_info,
                                          const SplitSolution& s, 
                                          SplitKKTResidual& kkt_residual) const {
  data.diff_3d = robot.CoM() - com_ref_;
  data.J_3d.setZero();
  robot.getCoMJacobian(data.J_3d);
  kkt_residual.lq().noalias() 
      += data.J_3d.transpose() * comf_weight_.asDiagonal() * data.diff_3d;
}


void CoMCost::evalTerminalCostHessian(Robot& robot, CostFunctionData& data, 
                                      const GridInfo& grid_info,
                                      const SplitSolution& s, 
                                      SplitKKTMatrix& kkt_matrix) const {
  data.J_3d.setZero();
  robot.getCoMJacobian(data.J_3d);
  kkt_matrix.Qqq().noalias()
      += data.J_3d.transpose() * comf_weight_.asDiagonal() * data.J_3d;
}


double CoMCost::evalImpulseCost(Robot& robot, const ImpulseStatus& impulse_status, 
                                CostFunctionData& data, 
                                const GridInfo& grid_info,
                                const ImpulseSplitSolution& s) const {
  double l = 0;
  data.diff_3d = robot.CoM() - com_ref_;
  l += (comi_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
  return 0.5 * l;
}


void CoMCost::evalImpulseCostDerivatives(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  data.J_3d.setZero();
  robot.getCoMJacobian(data.J_3d);
  kkt_residual.lq().noalias() 
      += data.J_3d.transpose() * comi_weight_.asDiagonal() * data.diff_3d;
}


void CoMCost::evalImpulseCostHessian(Robot& robot, 
                                     const ImpulseStatus& impulse_status, 
                                     CostFunctionData& data, 
                                     const GridInfo& grid_info,
                                     const ImpulseSplitSolution& s, 
                                     ImpulseSplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qqq().noalias()
      += data.J_3d.transpose() * comi_weight_.asDiagonal() * data.J_3d;
}

} // namespace robotoc