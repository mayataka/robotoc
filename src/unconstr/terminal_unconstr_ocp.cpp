#include "robotoc/unconstr/terminal_unconstr_ocp.hpp"

#include <stdexcept>
#include <cassert>


namespace robotoc {

TerminalUnconstrOCP::TerminalUnconstrOCP(
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


TerminalUnconstrOCP::TerminalUnconstrOCP() 
  : cost_(),
    constraints_(),
    contact_status_(),
    data_() {
}


bool TerminalUnconstrOCP::isFeasible(Robot& robot, const GridInfo& grid_info,
                                     const SplitSolution& s) {
  // return constraints_->isFeasible(robot, contact_status_, constraints_data_, s);
  return true;
}


void TerminalUnconstrOCP::initConstraints(Robot& robot, const GridInfo& grid_info, 
                                          const SplitSolution& s) { 
  // assert(time_step >= 0);
  // constraints_data_ = constraints_->createConstraintsData(robot, time_step);
  // constraints_->setSlackAndDual(robot, contact_status_, constraints_data_, s);
}


void TerminalUnconstrOCP::evalOCP(Robot& robot, const GridInfo& grid_info,
                                  const SplitSolution& s, 
                                  SplitKKTResidual& kkt_residual) {
  robot.updateKinematics(s.q);
  kkt_residual.setZero();
  data_.performance_index.cost = cost_->evalTerminalCost(robot, data_.cost_data, 
                                                         grid_info, s);
  // constraints_->evalConstraint(robot, contact_status_, data_.constraints_data, s);
  // data_.performance_index.cost_barrier = constraints_data_.logBarrier();
  // data_.performance_index.primal_feasibility 
  //     = data_.primalFeasibility<1>() + kkt_residual.primalFeasibility<1>();
  data_.performance_index.primal_feasibility = 0.0;
}


void TerminalUnconstrOCP::evalKKT(Robot& robot, const GridInfo& grid_info,
                                  const SplitSolution& s, 
                                  SplitKKTMatrix& kkt_matrix,
                                  SplitKKTResidual& kkt_residual) {
  robot.updateKinematics(s.q);
  kkt_matrix.setZero();
  kkt_residual.setZero();
  data_.performance_index.cost = cost_->quadratizeTerminalCost(robot, data_.cost_data, 
                                                               grid_info, s, 
                                                               kkt_residual, kkt_matrix);
  // constraints_->linearizeConstraints(robot, contact_status_, 
  //                                    data_.constraints_data, s, kkt_residual);
  // data_.performance_index.cost_barrier = constraints_data_.logBarrier();
  unconstr::stateequation::linearizeForwardEulerTerminal(s, kkt_residual);
  data_.performance_index.primal_feasibility = 0.0;
      // = data_.primalFeasibility<1>() + kkt_residual.primalFeasibility<1>();
  data_.performance_index.dual_feasibility
      = data_.dualFeasibility<1>() + kkt_residual.dualFeasibility<1>();
  data_.performance_index.kkt_error
      // = data_.KKTError() + kkt_residual.KKTError();
      = kkt_residual.KKTError();
  // constraints_->condenseSlackAndDual(contact_status_, constraints_data_, 
  //                                    kkt_matrix, kkt_residual);
}


void TerminalUnconstrOCP::expandPrimalAndDual(SplitDirection& d) {
  constraints_->expandSlackAndDual(contact_status_, data_.constraints_data, d);
}


double TerminalUnconstrOCP::maxPrimalStepSize() {
  return 1.0;
  // return constraints_->maxSlackStepSize(data_.constraints_data);
}


double TerminalUnconstrOCP::maxDualStepSize() {
  return 1.0;
  // return constraints_->maxDualStepSize(data_.constraints_data);
}


void TerminalUnconstrOCP::updatePrimal(const Robot& robot, 
                                       const double primal_step_size, 
                                       const SplitDirection& d, SplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  // constraints_->updateSlack(data_.constraints_data, primal_step_size);
}


void TerminalUnconstrOCP::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  // constraints_->updateDual(data_.constraints_data, dual_step_size);
}

} // namespace robotoc