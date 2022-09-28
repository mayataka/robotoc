#include "robotoc/ocp/terminal_stage.hpp"
#include "robotoc/dynamics/terminal_state_equation.hpp"

#include <cassert>


namespace robotoc {

TerminalStage::TerminalStage(const std::shared_ptr<CostFunction>& cost, 
                             const std::shared_ptr<Constraints>& constraints, 
                             const std::shared_ptr<ContactSequence>& contact_sequence)
  : cost_(cost), 
    constraints_(constraints),
    contact_sequence_(contact_sequence) {
}


TerminalStage::TerminalStage() 
  : cost_(), 
    constraints_(),
    contact_sequence_() {
}


OCPData TerminalStage::createData(const Robot& robot) const {
  OCPData data;
  data.performance_index = PerformanceIndex();
  data.cost_data = cost_->createCostFunctionData(robot);
  data.constraints_data = constraints_->createConstraintsData(robot);
  data.state_equation_data = StateEquationData(robot);
  data.contact_dynamics_data = ContactDynamicsData(robot);
  data.switching_constraint_data = SwitchingConstraintData(robot);
  return data;
}


bool TerminalStage::isFeasible(Robot& robot, const GridInfo& grid_info, 
                               const SplitSolution& s, OCPData& data) const {
  // TODO: add inequality constraints at the terminal OCP.
  return true;
}


void TerminalStage::initConstraints(Robot& robot, const GridInfo& grid_info, 
                                    const SplitSolution& s, OCPData& data) const {
  // TODO: add inequality constraints at the terminal OCP.
}


void TerminalStage::evalOCP(Robot& robot, const GridInfo& grid_info,
                            const SplitSolution& s, OCPData& data,
                            SplitKKTResidual& kkt_residual) const {
  // setup computation
  robot.updateKinematics(s.q, s.v);
  kkt_residual.setContactDimension(0);
  kkt_residual.setZero();
  // eval cost and constraints
  data.performance_index.cost = cost_->evalTerminalCost(robot, data.cost_data, grid_info, s);
  // constraints_->evalConstraint(robot, contact_status, data.constraints_data, s);
  // data.performance_index.cost_barrier = data.constraints_data.logBarrier();
  // summarize evaluations
  data.performance_index.primal_feasibility = data.primalFeasibility<1>();
}


void TerminalStage::evalKKT(Robot& robot, const GridInfo& grid_info,
                            const Eigen::VectorXd& q_prev, 
                            const SplitSolution& s, OCPData& data, 
                            SplitKKTMatrix& kkt_matrix, 
                            SplitKKTResidual& kkt_residual) const {
  // setup computation
  robot.updateKinematics(s.q, s.v);
  kkt_matrix.setContactDimension(0);
  kkt_residual.setContactDimension(0);
  kkt_matrix.setZero();
  kkt_residual.setZero();
  // eval cost and constraints
  data.performance_index.cost 
      = cost_->quadratizeTerminalCost(robot, data.cost_data, grid_info, s, 
                                      kkt_residual, kkt_matrix);
  // constraints_->linearizeConstraints(robot, contact_status, data.constraints_data, 
  //                                    s, kkt_residual);
  // data.performance_index.cost_barrier = data.constraints_data.logBarrier();
  // eval dynamics
  linearizeTerminalStateEquation(robot, q_prev, s, data.state_equation_data, 
                                 kkt_matrix, kkt_residual);
  // summarize evaluations
  data.performance_index.primal_feasibility = data.primalFeasibility<1>();
  data.performance_index.dual_feasibility 
      = data.dualFeasibility<1>() + kkt_residual.dualFeasibility<1>();
  data.performance_index.kkt_error = data.KKTError() + kkt_residual.KKTError();
  // Forms linear system
  // constraints_->condenseSlackAndDual(contact_status, data.constraints_data, 
  //                                    kkt_matrix, kkt_residual);
  correctLinearizeTerminalStateEquation(data.state_equation_data, kkt_matrix);
}


void TerminalStage::expandPrimal(const GridInfo& grid_info, OCPData& data, 
                                 SplitDirection& d) const {
  // const auto& contact_status = contact_sequence_->contactStatus(grid_info.contact_phase);
  // d.setContactDimension(contact_status.dimf());
  d.setContactDimension(0);
  // constraints_->expandSlackAndDual(contact_status, data.constraints_data, d);
}


void TerminalStage::expandDual(const GridInfo& grid_info, OCPData& data,
                               SplitDirection& d, const double dts) const {
  // assert(grid_info.dt > 0);
  correctCostateDirection(data.state_equation_data, d);
}


double TerminalStage::maxPrimalStepSize(const OCPData& data) const {
  // return constraints_->maxSlackStepSize(data.constraints_data);
  return 1.0;
}


double TerminalStage::maxDualStepSize(const OCPData& data) const {
  // return constraints_->maxDualStepSize(data.constraints_data);
  return 1.0;
}


void TerminalStage::updatePrimal(const Robot& robot, 
                                 const double primal_step_size, 
                                 const SplitDirection& d, SplitSolution& s, 
                                 OCPData& data) const {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  // constraints_->updateSlack(data.constraints_data, primal_step_size);
}


void TerminalStage::updateDual(const double dual_step_size, 
                               OCPData& data) const {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  // constraints_->updateDual(data.constraints_data, dual_step_size);
}

} // namespace robotoc