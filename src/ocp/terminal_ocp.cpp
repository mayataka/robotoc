#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/ocp/riccati_factorizer.hpp"

#include <cassert>


namespace idocp {

TerminalOCP::TerminalOCP(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(),
    kkt_residual_(robot),
    kkt_matrix_(robot),
    use_kinematics_(false) {
  if (cost_->useKinematics() || constraints_->useKinematics()) {
    use_kinematics_ = true;
  }
}


TerminalOCP::TerminalOCP() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    kkt_residual_(),
    kkt_matrix_(),
    use_kinematics_(false) {
}


TerminalOCP::~TerminalOCP() {
}


bool TerminalOCP::isFeasible(Robot& robot, const SplitSolution& s) {
  // TODO: add inequality constraints at the terminal OCP.
  return true;
}


void TerminalOCP::initConstraints(Robot& robot, const int time_step, 
                                  const double dtau, const SplitSolution& s) {
  assert(time_step >= 0);
  assert(dtau > 0);
  // TODO: add inequality constraints at the terminal OCP.
}


void TerminalOCP::linearizeOCP(Robot& robot, const double t, 
                               const SplitSolution& s) {
  kkt_residual_.lx().setZero();
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v);
  }
  cost_->computeTerminalCostDerivatives(robot, cost_data_, t, s, kkt_residual_);
  kkt_residual_.lq().noalias() -= s.lmd;
  kkt_residual_.lv().noalias() -= s.gmm;
  kkt_matrix_.Qqq().setZero();
  kkt_matrix_.Qvv().setZero();
  cost_->computeTerminalCostHessian(robot, cost_data_, t, s, kkt_matrix_);
}


void TerminalOCP::backwardRiccatiRecursion(
      RiccatiFactorization& riccati) const {
  riccati.sq = - kkt_residual_.lq();
  riccati.sv = - kkt_residual_.lv();
  riccati.Pqq = kkt_matrix_.Qqq();
  riccati.Pvv = kkt_matrix_.Qvv();
}


void TerminalOCP::computePrimalDirection(const RiccatiFactorization& riccati,
                                         SplitDirection& d) const {
  static constexpr bool exist_state_constraint = false;
  RiccatiFactorizer::computeCostateDirection(riccati, d, 
                                             exist_state_constraint);
}


void TerminalOCP::computeDualDirection(const SplitDirection& d) {
}

 
double TerminalOCP::maxPrimalStepSize() {
  return 1;
  // TODO: add inequality constraints at the terminal OCP.
}


double TerminalOCP::maxDualStepSize() {
  return 1;
  // TODO: add inequality constraints at the terminal OCP.
}


double TerminalOCP::terminalCost(Robot& robot, const double t, 
                                 const SplitSolution& s) {
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v);
  }
  return cost_->phi(robot, cost_data_, t, s);
}


void TerminalOCP::updateDual(const double step_size) {
  assert(step_size > 0);
  assert(step_size <= 1);
  // TODO: add inequality constraints at the terminal OCP.
}


void TerminalOCP::updatePrimal(Robot& robot, const double step_size, 
                               const SplitDirection& d,
                               SplitSolution& s) const {
  assert(step_size > 0);
  assert(step_size <= 1);
  s.lmd.noalias() += step_size * d.dlmd();
  s.gmm.noalias() += step_size * d.dgmm();
  robot.integrateConfiguration(d.dq(), step_size, s.q);
  s.v.noalias() += step_size * d.dv();
}


void TerminalOCP::computeKKTResidual(Robot& robot, const double t,  
                                     const SplitSolution& s) {

  kkt_residual_.lx().setZero();
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v);
  }
  cost_->computeTerminalCostDerivatives(robot, cost_data_, t, s, 
                                        kkt_residual_);
  kkt_residual_.lq().noalias() -= s.lmd;
  kkt_residual_.lv().noalias() -= s.gmm;
}


double TerminalOCP::squaredNormKKTResidual() const {
  double error = 0;
  error += kkt_residual_.lx().squaredNorm();
  return error;
}

} // namespace idocp