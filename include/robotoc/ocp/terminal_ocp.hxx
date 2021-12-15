#ifndef ROBOTOC_TERMINAL_OCP_HXX_
#define ROBOTOC_TERMINAL_OCP_HXX_

#include "robotoc/ocp/terminal_ocp.hpp"

#include <cassert>

namespace robotoc {

inline TerminalOCP::TerminalOCP(const Robot& robot, 
                                const std::shared_ptr<CostFunction>& cost, 
                                const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(),
    state_equation_(robot),
    use_kinematics_(false),
    terminal_cost_(0) {
  if (cost_->useKinematics() || constraints_->useKinematics()) {
    use_kinematics_ = true;
  }
}


inline TerminalOCP::TerminalOCP() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    state_equation_(),
    use_kinematics_(false),
    terminal_cost_(0) {
}


inline TerminalOCP::~TerminalOCP() {
}


inline bool TerminalOCP::isFeasible(Robot& robot, const SplitSolution& s) {
  // TODO: add inequality constraints at the terminal OCP.
  return true;
}


inline void TerminalOCP::initConstraints(Robot& robot, const int time_stage, 
                                         const SplitSolution& s) {
  assert(time_stage >= 0);
  // TODO: add inequality constraints at the terminal OCP.
}


inline void TerminalOCP::evalOCP(Robot& robot, const double t, 
                                 const SplitSolution& s, 
                                 SplitKKTResidual& kkt_residual) {
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v);
  }
  terminal_cost_ = cost_->evalTerminalCost(robot, cost_data_, t, s);
}


inline void TerminalOCP::computeKKTResidual(Robot& robot, const double t,  
                                            const Eigen::VectorXd& q_prev,
                                            const SplitSolution& s,
                                            SplitKKTMatrix& kkt_matrix,
                                            SplitKKTResidual& kkt_residual) {
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v);
  }
  kkt_residual.lx.setZero();
  terminal_cost_ = cost_->linearizeTerminalCost(robot, cost_data_, t, s, 
                                                kkt_residual);
  state_equation_.linearizeStateEquation(robot, q_prev, s, 
                                         kkt_matrix, kkt_residual);
  kkt_residual.kkt_error = KKTError(kkt_residual);
}


inline void TerminalOCP::computeKKTSystem(Robot& robot, const double t, 
                                          const Eigen::VectorXd& q_prev,
                                          const SplitSolution& s,
                                          SplitKKTMatrix& kkt_matrix, 
                                          SplitKKTResidual& kkt_residual) {
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v);
  }
  kkt_matrix.Qxx.setZero();
  kkt_residual.lx.setZero();
  terminal_cost_ = cost_->quadratizeTerminalCost(robot, cost_data_, t, s, 
                                                 kkt_residual, kkt_matrix);
  state_equation_.linearizeStateEquation(robot, q_prev, s, 
                                         kkt_matrix, kkt_residual);
  kkt_residual.kkt_error = KKTError(kkt_residual);
  state_equation_.correctLinearizedStateEquation(kkt_matrix);
}

 
inline double TerminalOCP::maxPrimalStepSize() {
  return 1;
  // TODO: add inequality constraints at the terminal OCP.
}


inline double TerminalOCP::maxDualStepSize() {
  return 1;
  // TODO: add inequality constraints at the terminal OCP.
}


inline void TerminalOCP::expandPrimal(const SplitSolution& s, 
                                      SplitDirection& d) {
  // TODO: add inequality constraints at the terminal OCP.
}


inline void TerminalOCP::expandDual(SplitDirection& d) {
  // TODO: add inequality constraints at the terminal OCP.
  state_equation_.correctCostateDirection(d);
}


inline void TerminalOCP::updatePrimal(const Robot& robot, 
                                      const double step_size, 
                                      const SplitDirection& d, 
                                      SplitSolution& s) const {
  assert(step_size > 0);
  assert(step_size <= 1);
  s.lmd.noalias() += step_size * d.dlmd();
  s.gmm.noalias() += step_size * d.dgmm();
  robot.integrateConfiguration(d.dq(), step_size, s.q);
  s.v.noalias() += step_size * d.dv();
}


inline void TerminalOCP::updateDual(const double step_size) {
  assert(step_size > 0);
  assert(step_size <= 1);
  // TODO: add inequality constraints at the terminal OCP.
}


inline double TerminalOCP::KKTError(const SplitKKTResidual& kkt_residual) {
  double err = 0;
  err += kkt_residual.lx.squaredNorm();
  return err;
}


inline double TerminalOCP::terminalCost() const {
  return terminal_cost_;
}

} // namespace robotoc

#endif // ROBOTOC_TERMINAL_OCP_HXX_ 