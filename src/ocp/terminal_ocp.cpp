#include "idocp/ocp/terminal_ocp.hpp"

#include <assert.h>


namespace idocp {

TerminalOCP::TerminalOCP(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(robot),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot)),
    kkt_residual_(robot),
    kkt_matrix_(robot),
    s_tmp_(robot),
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
    s_tmp_(),
    use_kinematics_(false) {
}


TerminalOCP::~TerminalOCP() {
}


bool TerminalOCP::isFeasible(const Robot& robot, const SplitSolution& s) {
  // TODO: add inequality constraints at the terminal OCP.
  return true;
}


void TerminalOCP::initConstraints(const Robot& robot, const int time_step, 
                                  const double dtau, const SplitSolution& s) {
  assert(time_step >= 0);
  assert(dtau > 0);
  // TODO: add inequality constraints at the terminal OCP.
}


void TerminalOCP::linearizeOCP(Robot& robot, const double t, 
                               const SplitSolution& s, 
                               RiccatiFactorization& riccati) {
  kkt_residual_.lq().setZero();
  kkt_residual_.lv().setZero();
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  cost_->computeTerminalCostDerivatives(robot, cost_data_, t, s, 
                                        kkt_residual_);
  riccati.sq = - kkt_residual_.lq() + s.lmd;
  riccati.sv = - kkt_residual_.lv() + s.gmm;
  kkt_matrix_.Qqq().setZero();
  kkt_matrix_.Qvv().setZero();
  cost_->computeTerminalCostHessian(robot, cost_data_, t, s, kkt_matrix_);
  riccati.Pqq = kkt_matrix_.Qqq();
  riccati.Pvv = kkt_matrix_.Qvv();
}


void TerminalOCP::computeCondensedDirection(Robot& robot, const double dtau, 
                                            SplitDirection& d) {
  // TODO: add inequality constraints at the terminal OCP.
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
  return cost_->phi(robot, cost_data_, t, s);
}


double TerminalOCP::terminalCost(Robot& robot, const double step_size, 
                                 const double t, const SplitSolution& s, 
                                 const SplitDirection& d) {
  assert(step_size > 0);
  assert(step_size <= 1);
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp_.q);
  s_tmp_.v = s.v + step_size * d.dv();
  return cost_->phi(robot, cost_data_, t, s_tmp_);
}


void TerminalOCP::updateDual(const double step_size) {
  assert(step_size > 0);
  assert(step_size <= 1);
  // TODO: add inequality constraints at the terminal OCP.
}


void TerminalOCP::updatePrimal(Robot& robot, const double step_size, 
                               const RiccatiFactorization& riccati,
                               const SplitDirection& d,
                               SplitSolution& s) const {
  assert(step_size > 0);
  assert(step_size <= 1);
  s.lmd.noalias() 
      += step_size * (riccati.Pqq * d.dq() + riccati.Pqv * d.dv() - riccati.sq);
  s.gmm.noalias() 
      += step_size * (riccati.Pvq * d.dq() + riccati.Pvv * d.dv() - riccati.sv);
  robot.integrateConfiguration(d.dq(), step_size, s.q);
  s.v.noalias() += step_size * d.dv();
}


double TerminalOCP::squaredKKTErrorNorm(Robot& robot, const double t, 
                                        const SplitSolution& s) const {
  double error = 0;
  error += (kkt_residual_.lq()-s.lmd).squaredNorm();
  error += (kkt_residual_.lv()-s.gmm).squaredNorm();
  return error;
}


double TerminalOCP::computeSquaredKKTErrorNorm(Robot& robot, const double t, 
                                               const SplitSolution& s) {
  kkt_residual_.lq().setZero();
  kkt_residual_.lv().setZero();
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  cost_->computeTerminalCostDerivatives(robot, cost_data_, t, s, 
                                        kkt_residual_);
  kkt_residual_.lq().noalias() -= s.lmd;
  kkt_residual_.lv().noalias() -= s.gmm;
  double error = 0;
  error += kkt_residual_.lq().squaredNorm();
  error += kkt_residual_.lv().squaredNorm();
  return error;
}

} // namespace idocp