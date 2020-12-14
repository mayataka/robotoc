#ifndef IDOCP_TERMINAL_OCP_HXX_
#define IDOCP_TERMINAL_OCP_HXX_

#include "idocp/ocp/terminal_ocp.hpp"

#include <cassert>

namespace idocp {

inline TerminalOCP::TerminalOCP(const Robot& robot, 
                                const std::shared_ptr<CostFunction>& cost, 
                                const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(),
    use_kinematics_(false) {
  if (cost_->useKinematics() || constraints_->useKinematics()) {
    use_kinematics_ = true;
  }
}


inline TerminalOCP::TerminalOCP() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    use_kinematics_(false) {
}


inline TerminalOCP::~TerminalOCP() {
}


inline bool TerminalOCP::isFeasible(Robot& robot, const SplitSolution& s) {
  // TODO: add inequality constraints at the terminal OCP.
  return true;
}


inline void TerminalOCP::initConstraints(Robot& robot, const int time_step, 
                                         const double dtau, 
                                         const SplitSolution& s) {
  assert(time_step >= 0);
  assert(dtau > 0);
  // TODO: add inequality constraints at the terminal OCP.
}


inline void TerminalOCP::linearizeOCP(Robot& robot, const double t, 
                                      const SplitSolution& s,
                                      SplitKKTMatrix& kkt_matrix, 
                                      SplitKKTResidual& kkt_residual) {
  kkt_residual.lx().setZero();
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v);
  }
  cost_->computeTerminalCostDerivatives(robot, cost_data_, t, s, kkt_residual);
  kkt_residual.lq().noalias() -= s.lmd;
  kkt_residual.lv().noalias() -= s.gmm;
  kkt_matrix.Qqq().setZero();
  kkt_matrix.Qvv().setZero();
  cost_->computeTerminalCostHessian(robot, cost_data_, t, s, kkt_matrix);
}

 
inline double TerminalOCP::maxPrimalStepSize() {
  return 1;
  // TODO: add inequality constraints at the terminal OCP.
}


inline double TerminalOCP::maxDualStepSize() {
  return 1;
  // TODO: add inequality constraints at the terminal OCP.
}


inline double TerminalOCP::terminalCost(Robot& robot, const double t, 
                                        const SplitSolution& s) {
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v);
  }
  return cost_->phi(robot, cost_data_, t, s);
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


inline void TerminalOCP::computeKKTResidual(Robot& robot, const double t,  
                                            const SplitSolution& s,
                                            SplitKKTResidual& kkt_residual) {

  kkt_residual.lx().setZero();
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v);
  }
  cost_->computeTerminalCostDerivatives(robot, cost_data_, t, s, 
                                        kkt_residual);
  kkt_residual.lq().noalias() -= s.lmd;
  kkt_residual.lv().noalias() -= s.gmm;
}


inline double TerminalOCP::squaredNormKKTResidual(
    const SplitKKTResidual& kkt_residual) const {
  double error = 0;
  error += kkt_residual.lx().squaredNorm();
  return error;
}

} // namespace idocp

#endif // IDOCP_TERMINAL_OCP_HXX_ 