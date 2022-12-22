#include "robotoc/unconstr/parnmpc_intermediate_stage.hpp"
#include "robotoc/dynamics/unconstr_state_equation.hpp"
#include "robotoc/dynamics/unconstr_dynamics.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

ParNMPCIntermediateStage::ParNMPCIntermediateStage(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    constraints_(constraints),
    contact_status_(robot.createContactStatus()) {
}


ParNMPCIntermediateStage::ParNMPCIntermediateStage() 
  : cost_(),
    constraints_(),
    contact_status_() {
}


UnconstrOCPData ParNMPCIntermediateStage::createData(const Robot& robot) const {
  UnconstrOCPData data;
  data.cost_data = cost_->createCostFunctionData(robot);
  data.constraints_data = constraints_->createConstraintsData(robot);
  data.unconstr_dynamics = UnconstrDynamics(robot);
  return data;
}


bool ParNMPCIntermediateStage::isFeasible(Robot& robot, 
                                          const GridInfo& grid_info,
                                          const SplitSolution& s, 
                                          UnconstrOCPData& data) const {
  return constraints_->isFeasible(robot, contact_status_, grid_info, s,
                                  data.constraints_data);
}


void ParNMPCIntermediateStage::initConstraints(Robot& robot, 
                                               const GridInfo& grid_info, 
                                               const SplitSolution& s,
                                               UnconstrOCPData& data) const { 
  data.constraints_data = constraints_->createConstraintsData(robot, grid_info.stage+1);
  constraints_->setSlackAndDual(robot, contact_status_, grid_info, s, data.constraints_data);
}



void ParNMPCIntermediateStage::evalOCP(Robot& robot, const GridInfo& grid_info, 
                                       const Eigen::VectorXd& q_prev, 
                                       const Eigen::VectorXd& v_prev, 
                                       const SplitSolution& s, 
                                       UnconstrOCPData& data,
                                       SplitKKTResidual& kkt_residual) const {
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  robot.updateKinematics(s.q);
  data.performance_index.setZero();
  kkt_residual.setZero();
  data.performance_index.cost = cost_->evalStageCost(robot, contact_status_, 
                                                     grid_info, s, data.cost_data);
  constraints_->evalConstraint(robot, contact_status_, grid_info, s,
                               data.constraints_data);
  data.performance_index.cost_barrier = data.constraints_data.logBarrier();
  evalUnconstrBackwardEuler(grid_info.dt, q_prev, v_prev, s, kkt_residual);
  data.unconstr_dynamics.evalUnconstrDynamics(robot, s);
  data.performance_index.primal_feasibility 
      = data.primalFeasibility<1>() + kkt_residual.primalFeasibility<1>();
}


void ParNMPCIntermediateStage::evalKKT(Robot& robot, const GridInfo& grid_info,
                                       const Eigen::VectorXd& q_prev, 
                                       const Eigen::VectorXd& v_prev, 
                                       const SplitSolution& s, 
                                       const SplitSolution& s_next,
                                       UnconstrOCPData& data, 
                                       SplitKKTMatrix& kkt_matrix,
                                       SplitKKTResidual& kkt_residual) const {
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  robot.updateKinematics(s.q);
  data.performance_index.setZero();
  kkt_matrix.setZero();
  kkt_residual.setZero();
  data.performance_index.cost = cost_->quadratizeStageCost(robot, contact_status_, grid_info,
                                                           s, data.cost_data, kkt_residual, kkt_matrix);
  constraints_->linearizeConstraints(robot, contact_status_, grid_info, s,
                                     data.constraints_data, kkt_residual);
  data.performance_index.cost_barrier = data.constraints_data.logBarrier();
  linearizeUnconstrBackwardEuler(grid_info.dt, q_prev, v_prev, s, s_next, 
                                 kkt_matrix, kkt_residual);
  data.unconstr_dynamics.linearizeUnconstrDynamics(robot, grid_info.dt, s, kkt_residual);
  data.performance_index.primal_feasibility 
      = data.primalFeasibility<1>() + kkt_residual.primalFeasibility<1>();
  data.performance_index.dual_feasibility
      = data.dualFeasibility<1>() + kkt_residual.dualFeasibility<1>();
  data.performance_index.kkt_error
      = data.KKTError() + kkt_residual.KKTError();
  constraints_->condenseSlackAndDual(contact_status_, grid_info,
                                     data.constraints_data, 
                                     kkt_matrix, kkt_residual);
  data.unconstr_dynamics.condenseUnconstrDynamics(kkt_matrix, kkt_residual);
}


void ParNMPCIntermediateStage::expandPrimalAndDual(
    const GridInfo& grid_info, const SplitKKTMatrix& kkt_matrix, 
    const SplitKKTResidual& kkt_residual, UnconstrOCPData& data, 
    SplitDirection& d) const {
  data.unconstr_dynamics.expandPrimal(d);
  data.unconstr_dynamics.expandDual(grid_info.dt, kkt_matrix, kkt_residual, d);
  constraints_->expandSlackAndDual(contact_status_, grid_info, d,
                                   data.constraints_data);
}


double ParNMPCIntermediateStage::maxPrimalStepSize(
    const UnconstrOCPData& data) const {
  return constraints_->maxSlackStepSize(data.constraints_data);
}


double ParNMPCIntermediateStage::maxDualStepSize(
    const UnconstrOCPData& data) const {
  return constraints_->maxDualStepSize(data.constraints_data);
}


void ParNMPCIntermediateStage::updatePrimal(const Robot& robot, 
                                            const double primal_step_size, 
                                            const SplitDirection& d, 
                                            SplitSolution& s,
                                            UnconstrOCPData& data) const {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  constraints_->updateSlack(data.constraints_data, primal_step_size);
}


void ParNMPCIntermediateStage::updateDual(const double dual_step_size,
                                          UnconstrOCPData& data) const {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(data.constraints_data, dual_step_size);
}

} // namespace robotoc