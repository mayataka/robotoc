#include "robotoc/ocp/terminal_ocp.hpp"
#include "robotoc/dynamics/terminal_state_equation.hpp"

#include <cassert>


namespace robotoc {

TerminalOCP::TerminalOCP(const Robot& robot, 
                         const std::shared_ptr<CostFunction>& cost, 
                         const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(),
    state_equation_data_(robot),
    terminal_cost_(0),
    barrier_cost_(0) {
}


TerminalOCP::TerminalOCP() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    state_equation_data_(),
    terminal_cost_(0),
    barrier_cost_(0) {
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
  terminal_cost_ = cost_->evalTerminalCost(robot, cost_data_, grid_info, s);
}


void TerminalOCP::computeKKTResidual(Robot& robot, const GridInfo& grid_info,
                                     const Eigen::VectorXd& q_prev,
                                     const SplitSolution& s,
                                     SplitKKTMatrix& kkt_matrix,
                                     SplitKKTResidual& kkt_residual) {
  robot.updateKinematics(s.q, s.v);
  kkt_residual.lx.setZero();
  terminal_cost_ = cost_->linearizeTerminalCost(robot, cost_data_, grid_info, s, 
                                                kkt_residual);
  TerminalStateEquation::linearize(robot, state_equation_data_, q_prev, s, 
                                   kkt_matrix, kkt_residual);
  kkt_residual.kkt_error = KKTError(kkt_residual);
}


void TerminalOCP::computeKKTSystem(Robot& robot, const GridInfo& grid_info,
                                   const Eigen::VectorXd& q_prev,
                                   const SplitSolution& s,
                                   SplitKKTMatrix& kkt_matrix, 
                                   SplitKKTResidual& kkt_residual) {
  robot.updateKinematics(s.q, s.v);
  kkt_matrix.Qxx.setZero();
  kkt_residual.lx.setZero();
  terminal_cost_ = cost_->quadratizeTerminalCost(robot, cost_data_, grid_info, s, 
                                                 kkt_residual, kkt_matrix);
  TerminalStateEquation::linearize(robot, state_equation_data_, q_prev, s, 
                                   kkt_matrix, kkt_residual);
  kkt_residual.kkt_error = KKTError(kkt_residual);
  TerminalStateEquation::correctLinearize(state_equation_data_, kkt_matrix);
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
  TerminalStateEquation::correctCostateDirection(state_equation_data_, d);
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
  return terminal_cost_;
}

} // namespace robotoc