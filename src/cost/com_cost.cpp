#include "idocp/cost/com_cost.hpp"


namespace idocp {

CoMCost::CoMCost(const Robot& robot)
  : CostFunctionComponentBase(),
    CoM_ref_(Eigen::Vector3d::Zero()),
    q_weight_(Eigen::Vector3d::Zero()),
    qf_weight_(Eigen::Vector3d::Zero()),
    qi_weight_(Eigen::Vector3d::Zero()) {
}


CoMCost::CoMCost()
  : CostFunctionComponentBase(),
    CoM_ref_(),
    q_weight_(),
    qf_weight_(),
    qi_weight_() {
}


CoMCost::~CoMCost() {
}


void CoMCost::set_CoM_ref(const Eigen::Vector3d& CoM_ref) {
  CoM_ref_ = CoM_ref;
}


void CoMCost::set_q_weight(const Eigen::Vector3d& q_weight) {
  q_weight_ = q_weight;
}


void CoMCost::set_qf_weight(const Eigen::Vector3d& qf_weight) {
  qf_weight_ = qf_weight;
}


void CoMCost::set_qi_weight(const Eigen::Vector3d& qi_weight) {
  qi_weight_ = qi_weight;
}


bool CoMCost::useKinematics() const {
  return true;
}


double CoMCost::computeStageCost(Robot& robot, CostFunctionData& data, 
                                 const double t, const double dt, 
                                 const SplitSolution& s) const {
  double l = 0;
  data.diff_3d = robot.CoM() - CoM_ref_;
  l += (q_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
  return 0.5 * dt * l;
}


void CoMCost::computeStageCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, const double dt, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  data.J_3d.setZero();
  robot.getCoMJacobian(data.J_3d);
  kkt_residual.lq().noalias() 
      += dt * data.J_3d.transpose() * q_weight_.asDiagonal() * data.diff_3d;
}


void CoMCost::computeStageCostHessian(
    Robot& robot, CostFunctionData& data, const double t, const double dt, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qqq().noalias()
      += dt * data.J_3d.transpose() * q_weight_.asDiagonal() * data.J_3d;
}


double CoMCost::computeTerminalCost(Robot& robot, CostFunctionData& data, 
                                    const double t, 
                                    const SplitSolution& s) const {
  double l = 0;
  data.diff_3d = robot.CoM() - CoM_ref_;
  l += (qf_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
  return 0.5 * l;
}


void CoMCost::computeTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  data.diff_3d = robot.CoM() - CoM_ref_;
  data.J_3d.setZero();
  robot.getCoMJacobian(data.J_3d);
  kkt_residual.lq().noalias() 
      += data.J_3d.transpose() * qf_weight_.asDiagonal() * data.diff_3d;
}


void CoMCost::computeTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  data.J_3d.setZero();
  robot.getCoMJacobian(data.J_3d);
  kkt_matrix.Qqq().noalias()
      += data.J_3d.transpose() * qf_weight_.asDiagonal() * data.J_3d;
}


double CoMCost::computeImpulseCost(Robot& robot,  
                                           CostFunctionData& data, 
                                           const double t, 
                                           const ImpulseSplitSolution& s) const {
  double l = 0;
  data.diff_3d = robot.CoM() - CoM_ref_;
  l += (qi_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
  return 0.5 * l;
}


void CoMCost::computeImpulseCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  data.J_3d.setZero();
  robot.getCoMJacobian(data.J_3d);
  kkt_residual.lq().noalias() 
      += data.J_3d.transpose() * qi_weight_.asDiagonal() * data.diff_3d;
}


void CoMCost::computeImpulseCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qqq().noalias()
      += data.J_3d.transpose() * qi_weight_.asDiagonal() * data.J_3d;
}

} // namespace idocp