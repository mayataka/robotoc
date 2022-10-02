#include "robotoc/unconstr/unconstr_intermediate_stage.hpp"

#include <stdexcept>
#include <cassert>


namespace robotoc {

UnconstrIntermediateStage::UnconstrIntermediateStage(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    constraints_(constraints),
    contact_status_(robot.createContactStatus()) {
}


UnconstrIntermediateStage::UnconstrIntermediateStage() 
  : cost_(),
    constraints_(),
    contact_status_() {
}


UnconstrOCPData UnconstrIntermediateStage::createData(const Robot& robot) const {
  UnconstrOCPData data;
  data.cost_data = cost_->createCostFunctionData(robot);
  data.constraints_data = constraints_->createConstraintsData(robot);
  data.unconstr_dynamics = UnconstrDynamics(robot);
  return data;
}


bool UnconstrIntermediateStage::isFeasible(Robot& robot, 
                                           const GridInfo& grid_info,
                                           const SplitSolution& s,
                                           UnconstrOCPData& data) const {
  return constraints_->isFeasible(robot, contact_status_, data.constraints_data, s);
}


void UnconstrIntermediateStage::initConstraints(Robot& robot, 
                                                const GridInfo& grid_info, 
                                                const SplitSolution& s,
                                                UnconstrOCPData& data) const {
  data.constraints_data = constraints_->createConstraintsData(robot, grid_info.time_stage);
  constraints_->setSlackAndDual(robot, contact_status_, data.constraints_data, s);
}


void UnconstrIntermediateStage::evalOCP(Robot& robot, const GridInfo& grid_info,
                                        const SplitSolution& s, 
                                        const SplitSolution& s_next, 
                                        UnconstrOCPData& data, 
                                        SplitKKTResidual& kkt_residual) const {
  robot.updateKinematics(s.q);
  kkt_residual.setZero();
  data.performance_index.cost = cost_->evalStageCost(robot, contact_status_, 
                                                      data.cost_data, grid_info, s);
  constraints_->evalConstraint(robot, contact_status_, data.constraints_data, s);
  data.performance_index.cost_barrier = data.constraints_data.logBarrier();
  unconstr::stateequation::evalForwardEuler(grid_info.dt, s, s_next, kkt_residual);
  data.unconstr_dynamics.evalUnconstrDynamics(robot, s);
  data.performance_index.primal_feasibility 
      = data.primalFeasibility<1>() + kkt_residual.primalFeasibility<1>();
}


void UnconstrIntermediateStage::evalKKT(Robot& robot, const GridInfo& grid_info,
                                        const SplitSolution& s, 
                                        const SplitSolution& s_next, 
                                        UnconstrOCPData& data, 
                                        SplitKKTMatrix& kkt_matrix,
                                        SplitKKTResidual& kkt_residual) const {
  robot.updateKinematics(s.q);
  kkt_matrix.setZero();
  kkt_residual.setZero();
  data.performance_index.cost = cost_->quadratizeStageCost(robot, contact_status_, 
                                                           data.cost_data, grid_info, s, 
                                                           kkt_residual, kkt_matrix);
  constraints_->linearizeConstraints(robot, contact_status_, 
                                     data.constraints_data, s, kkt_residual);
  data.performance_index.cost_barrier = data.constraints_data.logBarrier();
  unconstr::stateequation::linearizeForwardEuler(grid_info.dt, s, s_next, 
                                                 kkt_matrix, kkt_residual);
  data.unconstr_dynamics.linearizeUnconstrDynamics(robot, grid_info.dt, s, kkt_residual);
  data.performance_index.primal_feasibility 
      = data.primalFeasibility<1>() + kkt_residual.primalFeasibility<1>();
  data.performance_index.dual_feasibility
      = data.dualFeasibility<1>() + kkt_residual.dualFeasibility<1>();
  data.performance_index.kkt_error
      = data.KKTError() + kkt_residual.KKTError();
  constraints_->condenseSlackAndDual(contact_status_, data.constraints_data, 
                                     kkt_matrix, kkt_residual);
  data.unconstr_dynamics.condenseUnconstrDynamics(kkt_matrix, kkt_residual);
}


void UnconstrIntermediateStage::expandPrimalAndDual(
    const double dt, const SplitKKTMatrix& kkt_matrix, 
    const SplitKKTResidual& kkt_residual, UnconstrOCPData& data, 
    SplitDirection& d) const {
  assert(dt > 0);
  data.unconstr_dynamics.expandPrimal(d);
  data.unconstr_dynamics.expandDual(dt, kkt_matrix, kkt_residual, d);
  constraints_->expandSlackAndDual(contact_status_, data.constraints_data, d);
}


double UnconstrIntermediateStage::maxPrimalStepSize(
    const UnconstrOCPData& data) const {
  return constraints_->maxSlackStepSize(data.constraints_data);
}


double UnconstrIntermediateStage::maxDualStepSize(
    const UnconstrOCPData& data) const {
  return constraints_->maxDualStepSize(data.constraints_data);
}


void UnconstrIntermediateStage::updatePrimal(const Robot& robot, 
                                             const double primal_step_size, 
                                             const SplitDirection& d, 
                                             SplitSolution& s,
                                             UnconstrOCPData& data) const {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  constraints_->updateSlack(data.constraints_data, primal_step_size);
}


void UnconstrIntermediateStage::updateDual(const double dual_step_size,
                                           UnconstrOCPData& data) const {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(data.constraints_data, dual_step_size);
}

} // namespace robotoc