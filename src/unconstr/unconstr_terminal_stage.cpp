#include "robotoc/unconstr/unconstr_terminal_stage.hpp"

#include <stdexcept>
#include <cassert>


namespace robotoc {

UnconstrTerminalStage::UnconstrTerminalStage(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    constraints_(constraints),
    contact_status_(robot.createContactStatus()) {
}


UnconstrTerminalStage::UnconstrTerminalStage() 
  : cost_(),
    constraints_(),
    contact_status_() {
}


UnconstrOCPData UnconstrTerminalStage::createData(const Robot& robot) const {
  UnconstrOCPData data;
  data.cost_data = cost_->createCostFunctionData(robot);
  data.constraints_data = constraints_->createConstraintsData(robot);
  data.unconstr_dynamics = UnconstrDynamics(robot);
  return data;
}


bool UnconstrTerminalStage::isFeasible(Robot& robot, const GridInfo& grid_info,
                                       const SplitSolution& s, 
                                       UnconstrOCPData& data) const {
  // return constraints_->isFeasible(robot, contact_status_, constraints_data_, s);
  return true;
}


void UnconstrTerminalStage::initConstraints(Robot& robot, const GridInfo& grid_info, 
                                            const SplitSolution& s, 
                                            UnconstrOCPData& data) const { 
  // assert(time_step >= 0);
  // constraints_data_ = constraints_->createConstraintsData(robot, time_step);
  // constraints_->setSlackAndDual(robot, contact_status_, constraints_data_, s);
}


void UnconstrTerminalStage::evalOCP(Robot& robot, const GridInfo& grid_info,
                                    const SplitSolution& s, 
                                    UnconstrOCPData& data, 
                                    SplitKKTResidual& kkt_residual) const {
  robot.updateKinematics(s.q);
  kkt_residual.setZero();
  data.performance_index.cost = cost_->evalTerminalCost(robot, data.cost_data, 
                                                         grid_info, s);
  // constraints_->evalConstraint(robot, contact_status_, data.constraints_data, s);
  // data.performance_index.cost_barrier = constraints_data_.logBarrier();
  // data.performance_index.primal_feasibility 
  //     = data.primalFeasibility<1>() + kkt_residual.primalFeasibility<1>();
  data.performance_index.primal_feasibility = 0.0;
}


void UnconstrTerminalStage::evalKKT(Robot& robot, const GridInfo& grid_info,
                                    const SplitSolution& s, 
                                    UnconstrOCPData& data, 
                                    SplitKKTMatrix& kkt_matrix,
                                    SplitKKTResidual& kkt_residual) const {
  robot.updateKinematics(s.q);
  kkt_matrix.setZero();
  kkt_residual.setZero();
  data.performance_index.cost = cost_->quadratizeTerminalCost(robot, data.cost_data, 
                                                               grid_info, s, 
                                                               kkt_residual, kkt_matrix);
  // constraints_->linearizeConstraints(robot, contact_status_, 
  //                                    data.constraints_data, s, kkt_residual);
  // data.performance_index.cost_barrier = constraints_data_.logBarrier();
  unconstr::stateequation::linearizeForwardEulerTerminal(s, kkt_residual);
  data.performance_index.primal_feasibility = 0.0;
      // = data.primalFeasibility<1>() + kkt_residual.primalFeasibility<1>();
  data.performance_index.dual_feasibility
      = data.dualFeasibility<1>() + kkt_residual.dualFeasibility<1>();
  data.performance_index.kkt_error
      // = data.KKTError() + kkt_residual.KKTError();
      = kkt_residual.KKTError();
  // constraints_->condenseSlackAndDual(contact_status_, constraints_data_, 
  //                                    kkt_matrix, kkt_residual);
}


void UnconstrTerminalStage::expandPrimalAndDual(UnconstrOCPData& data, 
                                                SplitDirection& d) const {
  // constraints_->expandSlackAndDual(contact_status_, data.constraints_data, d);
}


double UnconstrTerminalStage::maxPrimalStepSize(const UnconstrOCPData& data) const {
  return 1.0;
  // return constraints_->maxSlackStepSize(data.constraints_data);
}


double UnconstrTerminalStage::maxDualStepSize(const UnconstrOCPData& data) const {
  return 1.0;
  // return constraints_->maxDualStepSize(data.constraints_data);
}


void UnconstrTerminalStage::updatePrimal(const Robot& robot, 
                                         const double primal_step_size, 
                                         const SplitDirection& d, 
                                         SplitSolution& s,
                                         UnconstrOCPData& data) const {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  // constraints_->updateSlack(data.constraints_data, primal_step_size);
}


void UnconstrTerminalStage::updateDual(const double dual_step_size, 
                                       UnconstrOCPData& data) const {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  // constraints_->updateDual(data.constraints_data, dual_step_size);
}

} // namespace robotoc