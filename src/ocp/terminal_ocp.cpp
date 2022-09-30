#include "robotoc/ocp/terminal_ocp.hpp"
#include "robotoc/dynamics/terminal_state_equation.hpp"

#include <cassert>


namespace robotoc {

TerminalOCP::TerminalOCP(const Robot& robot, 
                         const std::shared_ptr<CostFunction>& cost, 
                         const std::shared_ptr<Constraints>& constraints) {
  ocp_.cost = cost;
  ocp_.constraints = constraints;
  data_.cost_data = cost->createCostFunctionData(robot);
  data_.constraints_data = constraints->createConstraintsData(robot, 0);
  data_.state_equation_data = StateEquationData(robot);
  data_.contact_dynamics_data = ContactDynamicsData(robot);
  data_.switching_constraint_data = SwitchingConstraintData(robot);
}


TerminalOCP::TerminalOCP() 
  : ocp_(),
    data_() {
}


TerminalOCP::~TerminalOCP() {
}


bool TerminalOCP::isFeasible(Robot& robot, const SplitSolution& s) {
  // TODO: add inequality constraints at the terminal OCP.
  return true;
}


void TerminalOCP::initConstraints(Robot& robot, const int time_stage, 
                                  const SplitSolution& s) {
  assert(time_stage >= 0);
  // TODO: add inequality constraints at the terminal OCP.
}


void TerminalOCP::evalOCP(Robot& robot, const GridInfo& grid_info,
                          const SplitSolution& s, 
                          SplitKKTResidual& kkt_residual) {
  robot.updateKinematics(s.q, s.v);
  data_.performance_index.cost = ocp_.cost->evalTerminalCost(robot, data_.cost_data, grid_info, s);
}


void TerminalOCP::computeKKTResidual(Robot& robot, const GridInfo& grid_info,
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
}


void TerminalOCP::computeKKTSystem(Robot& robot, const GridInfo& grid_info,
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
  correctLinearizeTerminalStateEquation(data_.state_equation_data, kkt_matrix);
}

 
double TerminalOCP::maxPrimalStepSize() {
  return 1;
  // TODO: add inequality constraints at the terminal OCP.
}


double TerminalOCP::maxDualStepSize() {
  return 1;
  // TODO: add inequality constraints at the terminal OCP.
}


void TerminalOCP::expandPrimal(const SplitSolution& s, SplitDirection& d) {
  // TODO: add inequality constraints at the terminal OCP.
}


void TerminalOCP::expandDual(SplitDirection& d) {
  // TODO: add inequality constraints at the terminal OCP.
  correctCostateDirection(data_.state_equation_data, d);
}


void TerminalOCP::updatePrimal(const Robot& robot, const double step_size, 
                               const SplitDirection& d, 
                               SplitSolution& s) const {
  assert(step_size > 0);
  assert(step_size <= 1);
  s.lmd.noalias() += step_size * d.dlmd();
  s.gmm.noalias() += step_size * d.dgmm();
  robot.integrateConfiguration(d.dq(), step_size, s.q);
  s.v.noalias() += step_size * d.dv();
}


void TerminalOCP::updateDual(const double step_size) {
  assert(step_size > 0);
  assert(step_size <= 1);
  // TODO: add inequality constraints at the terminal OCP.
}


double TerminalOCP::KKTError(const SplitKKTResidual& kkt_residual) {
  double err = 0;
  err += kkt_residual.lx.squaredNorm();
  return err;
}


double TerminalOCP::terminalCost(const bool) const {
  return data_.performance_index.cost;
}

} // namespace robotoc