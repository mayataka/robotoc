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
                            const SplitSolution& s, 
                            SplitKKTResidual& kkt_residual) {
  robot.updateKinematics(s.q, s.v);
  data_.performance_index.cost = ocp_.cost->evalTerminalCost(robot, data_.cost_data, grid_info, s);
}


void TerminalStage::computeKKTResidual(Robot& robot, const GridInfo& grid_info,
                                     const Eigen::VectorXd& q_prev,
                                     const SplitSolution& s,
                                     SplitKKTMatrix& kkt_matrix,
                                     SplitKKTResidual& kkt_residual) {
  robot.updateKinematics(s.q, s.v);
  kkt_residual.lx.setZero();
  data_.performance_index.cost = ocp_.cost->linearizeTerminalCost(robot, data_.cost_data, grid_info, s, 
                                                kkt_residual);
  linearizeTerminalStateEquation(robot, q_prev, s, data_.state_equation_data, 
                                 kkt_matrix, kkt_residual);
  kkt_residual.kkt_error = KKTError(kkt_residual);
}


void TerminalStage::computeKKTSystem(Robot& robot, const GridInfo& grid_info,
                                   const Eigen::VectorXd& q_prev,
                                   const SplitSolution& s,
                                   SplitKKTMatrix& kkt_matrix, 
                                   SplitKKTResidual& kkt_residual) {
  robot.updateKinematics(s.q, s.v);
  kkt_matrix.Qxx.setZero();
  kkt_residual.lx.setZero();
  data_.performance_index.cost = ocp_.cost->quadratizeTerminalCost(robot, data_.cost_data, grid_info, s, 
                                                 kkt_residual, kkt_matrix);
  linearizeTerminalStateEquation(robot, q_prev, s, data_.state_equation_data, 
                                 kkt_matrix, kkt_residual);
  kkt_residual.kkt_error = KKTError(kkt_residual);
  correctLinearizeTerminalStateEquation(data_.state_equation_data, kkt_matrix);
}

 
double TerminalStage::maxPrimalStepSize() {
  return 1;
  // TODO: add inequality constraints at the terminal OCP.
}


double TerminalStage::maxDualStepSize() {
  return 1;
  // TODO: add inequality constraints at the terminal OCP.
}


void TerminalStage::expandPrimal(const SplitSolution& s, SplitDirection& d) {
  // TODO: add inequality constraints at the terminal OCP.
}


void TerminalStage::expandDual(SplitDirection& d) {
  // TODO: add inequality constraints at the terminal OCP.
  correctCostateDirection(data_.state_equation_data, d);
}


void TerminalStage::updatePrimal(const Robot& robot, const double step_size, 
                               const SplitDirection& d, 
                               SplitSolution& s) const {
  assert(step_size > 0);
  assert(step_size <= 1);
  s.lmd.noalias() += step_size * d.dlmd();
  s.gmm.noalias() += step_size * d.dgmm();
  robot.integrateConfiguration(d.dq(), step_size, s.q);
  s.v.noalias() += step_size * d.dv();
}


void TerminalStage::updateDual(const double step_size) {
  assert(step_size > 0);
  assert(step_size <= 1);
  // TODO: add inequality constraints at the terminal OCP.
}


double TerminalStage::KKTError(const SplitKKTResidual& kkt_residual) {
  double err = 0;
  err += kkt_residual.lx.squaredNorm();
  return err;
}


double TerminalStage::terminalCost(const bool) const {
  return data_.performance_index.cost;
}

} // namespace robotoc