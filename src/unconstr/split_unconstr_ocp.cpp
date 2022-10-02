#include "robotoc/unconstr/split_unconstr_ocp.hpp"

#include <stdexcept>
#include <cassert>


namespace robotoc {

SplitUnconstrOCP::SplitUnconstrOCP(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    constraints_(constraints),
    contact_status_(robot.createContactStatus()),
    data_() {
  data_.cost_data = cost->createCostFunctionData(robot);
  data_.constraints_data = constraints->createConstraintsData(robot);
  data_.unconstr_dynamics = UnconstrDynamics(robot);
}


SplitUnconstrOCP::SplitUnconstrOCP() 
  : cost_(),
    constraints_(),
    contact_status_(),
    data_() {
}


bool SplitUnconstrOCP::isFeasible(Robot& robot, const GridInfo& grid_info,
                                  const SplitSolution& s) {
  return constraints_->isFeasible(robot, contact_status_, data_.constraints_data, s);
}


void SplitUnconstrOCP::initConstraints(Robot& robot, const GridInfo& grid_info, 
                                       const SplitSolution& s) {
  data_.constraints_data = constraints_->createConstraintsData(robot, grid_info.time_stage);
  constraints_->setSlackAndDual(robot, contact_status_, data_.constraints_data, s);
}


void SplitUnconstrOCP::evalOCP(Robot& robot, const GridInfo& grid_info,
                               const SplitSolution& s, 
                               const SplitSolution& s_next, 
                               SplitKKTResidual& kkt_residual) {
  robot.updateKinematics(s.q);
  kkt_residual.setZero();
  data_.performance_index.cost = cost_->evalStageCost(robot, contact_status_, 
                                                      data_.cost_data, grid_info, s);
  constraints_->evalConstraint(robot, contact_status_, data_.constraints_data, s);
  data_.performance_index.cost_barrier = data_.constraints_data.logBarrier();
  unconstr::stateequation::evalForwardEuler(grid_info.dt, s, s_next, kkt_residual);
  data_.unconstr_dynamics.evalUnconstrDynamics(robot, s);
  data_.performance_index.primal_feasibility 
      = data_.primalFeasibility<1>() + kkt_residual.primalFeasibility<1>();
}


void SplitUnconstrOCP::evalKKT(Robot& robot, const GridInfo& grid_info,
                               const SplitSolution& s, 
                               const SplitSolution& s_next, 
                               SplitKKTMatrix& kkt_matrix,
                               SplitKKTResidual& kkt_residual) {
  robot.updateKinematics(s.q);
  kkt_matrix.setZero();
  kkt_residual.setZero();
  data_.performance_index.cost = cost_->quadratizeStageCost(robot, contact_status_, 
                                                            data_.cost_data, grid_info, s, 
                                                            kkt_residual, kkt_matrix);
  constraints_->linearizeConstraints(robot, contact_status_, 
                                     data_.constraints_data, s, kkt_residual);
  data_.performance_index.cost_barrier = data_.constraints_data.logBarrier();
  unconstr::stateequation::linearizeForwardEuler(grid_info.dt, s, s_next, 
                                                 kkt_matrix, kkt_residual);
  data_.unconstr_dynamics.linearizeUnconstrDynamics(robot, grid_info.dt, s, kkt_residual);
  data_.performance_index.primal_feasibility 
      = data_.primalFeasibility<1>() + kkt_residual.primalFeasibility<1>();
  data_.performance_index.dual_feasibility
      = data_.dualFeasibility<1>() + kkt_residual.dualFeasibility<1>();
  data_.performance_index.kkt_error
      = data_.KKTError() + kkt_residual.KKTError();
  constraints_->condenseSlackAndDual(contact_status_, data_.constraints_data, 
                                     kkt_matrix, kkt_residual);
  data_.unconstr_dynamics.condenseUnconstrDynamics(kkt_matrix, kkt_residual);
}


void SplitUnconstrOCP::expandPrimalAndDual(const double dt, 
                                           const SplitKKTMatrix& kkt_matrix, 
                                           const SplitKKTResidual& kkt_residual, 
                                           SplitDirection& d) {
  assert(dt > 0);
  data_.unconstr_dynamics.expandPrimal(d);
  data_.unconstr_dynamics.expandDual(dt, kkt_matrix, kkt_residual, d);
  constraints_->expandSlackAndDual(contact_status_, data_.constraints_data, d);
}


double SplitUnconstrOCP::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(data_.constraints_data);
}


double SplitUnconstrOCP::maxDualStepSize() {
  return constraints_->maxDualStepSize(data_.constraints_data);
}


void SplitUnconstrOCP::updatePrimal(const Robot& robot, 
                                    const double primal_step_size, 
                                    const SplitDirection& d, SplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  constraints_->updateSlack(data_.constraints_data, primal_step_size);
}


void SplitUnconstrOCP::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(data_.constraints_data, dual_step_size);
}

} // namespace robotoc