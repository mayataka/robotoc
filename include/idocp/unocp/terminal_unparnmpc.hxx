#ifndef IDOCP_TERMINAL_UNPARNMPC_HXX_ 
#define IDOCP_TERMINAL_UNPARNMPC_HXX_

#include "idocp/unocp/terminal_unparnmpc.hpp"

#include <stdexcept>
#include <cassert>

namespace idocp {

inline TerminalUnParNMPC::TerminalUnParNMPC(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot, 0)),
    unconstrained_dynamics_(robot),
    use_kinematics_(false),
    kkt_matrix_(robot),
    kkt_residual_(robot) {
  if (cost_->useKinematics() || constraints_->useKinematics() 
                             || robot.maxPointContacts() > 0) {
    use_kinematics_ = true;
  }
  try {
    if (robot.hasFloatingBase()) {
      throw std::logic_error("robot has floating base: robot should have no constraints!");
    }
    if (robot.maxPointContacts() > 0) {
      throw std::logic_error("robot can have contacts: robot should have no constraints!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline TerminalUnParNMPC::TerminalUnParNMPC() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    unconstrained_dynamics_(),
    use_kinematics_(false),
    kkt_matrix_(),
    kkt_residual_() {
}


inline TerminalUnParNMPC::~TerminalUnParNMPC() {
}


inline bool TerminalUnParNMPC::isFeasible(Robot& robot, const SplitSolution& s) {
  return constraints_->isFeasible(robot, constraints_data_, s);
}


inline void TerminalUnParNMPC::initConstraints(Robot& robot, 
                                               const int time_step, 
                                               const SplitSolution& s) { 
  assert(time_step >= 0);
  constraints_data_ = constraints_->createConstraintsData(robot, time_step);
  constraints_->setSlackAndDual(robot, constraints_data_, s);
}


inline void TerminalUnParNMPC::linearizeOCP(Robot& robot, const double t, 
                                            const double dtau, 
                                            const Eigen::VectorXd& q_prev, 
                                            const Eigen::VectorXd& v_prev, 
                                            const SplitSolution& s, 
                                            SplitUnKKTMatrix& unkkt_matrix,
                                            SplitUnKKTResidual& unkkt_residual) {
  assert(dtau >= 0);
  assert(q_prev.size() == robot.dimq());
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  kkt_matrix_.Qqq().setZero();
  kkt_matrix_.Qvv().diagonal().setZero();
  kkt_matrix_.Qaa().diagonal().setZero();
  kkt_matrix_.Quu().diagonal().setZero();
  kkt_residual_.setZero();
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  cost_->computeTerminalCostDerivatives(robot, cost_data_, t, s, kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s,
                                    kkt_residual_);
  stateequation::LinearizeBackwardEulerTerminal(robot, dtau, q_prev, v_prev, s, 
                                                kkt_matrix_, kkt_residual_);
  unconstrained_dynamics_.linearizeUnconstrainedDynamics(robot, dtau, s, 
                                                         kkt_residual_);
  cost_->computeStageCostHessian(robot, cost_data_, t, dtau, s, kkt_matrix_);
  cost_->computeTerminalCostHessian(robot, cost_data_, t, s, kkt_matrix_);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dtau, s, 
                                     kkt_matrix_, kkt_residual_);
  unconstrained_dynamics_.condenseUnconstrainedDynamics(
      kkt_matrix_, kkt_residual_, unkkt_matrix, unkkt_residual);
}


inline void TerminalUnParNMPC::computeCondensedDirection(Robot& robot, 
                                                         const double dtau, 
                                                         const SplitSolution& s, 
                                                         SplitDirection& d) {
  unconstrained_dynamics_.computeCondensedDirection(dtau, kkt_matrix_, 
                                                    kkt_residual_, d);
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, s, d);
}


inline double TerminalUnParNMPC::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


inline double TerminalUnParNMPC::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


inline void TerminalUnParNMPC::updatePrimal(const Robot& robot, 
                                            const double primal_step_size, 
                                            const SplitDirection& d, 
                                            SplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  constraints_->updateSlack(constraints_data_, primal_step_size);
}


inline void TerminalUnParNMPC::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(constraints_data_, dual_step_size);
}


inline void TerminalUnParNMPC::computeKKTResidual(Robot& robot, const double t, 
                                                  const double dtau, 
                                                  const Eigen::VectorXd& q_prev, 
                                                  const Eigen::VectorXd& v_prev, 
                                                  const SplitSolution& s) {
  assert(dtau >= 0);
  assert(q_prev.size() == robot.dimq());
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  kkt_residual_.setZero();
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  cost_->computeTerminalCostDerivatives(robot, cost_data_, t, s, kkt_residual_);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s,
                                    kkt_residual_);
  stateequation::LinearizeBackwardEulerTerminal(robot, dtau, q_prev, v_prev, s, 
                                                kkt_matrix_, kkt_residual_);
  unconstrained_dynamics_.linearizeUnconstrainedDynamics(robot, dtau, s, 
                                                         kkt_residual_);
}


inline double TerminalUnParNMPC::squaredNormKKTResidual(
    const double dtau) const {
  double error = 0;
  error += kkt_residual_.lx().squaredNorm();
  error += kkt_residual_.la.squaredNorm();
  error += kkt_residual_.lu().squaredNorm();
  error += stateequation::SquaredNormStateEuqationResidual(kkt_residual_);
  error += unconstrained_dynamics_.squaredNormUnconstrainedDynamicsResidual(dtau);
  error += dtau * dtau * constraints_->squaredNormPrimalAndDualResidual(constraints_data_);
  return error;
}


inline double TerminalUnParNMPC::stageCost(Robot& robot, const double t,  
                                           const double dtau, 
                                           const SplitSolution& s, 
                                           const double primal_step_size) {
  assert(dtau >= 0);
  assert(primal_step_size >= 0);
  assert(primal_step_size <= 1);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  double cost = 0;
  cost += cost_->computeStageCost(robot, cost_data_, t, dtau, s);
  cost += cost_->computeTerminalCost(robot, cost_data_, t, s);
  if (primal_step_size > 0) {
    cost += dtau * constraints_->costSlackBarrier(constraints_data_, 
                                                  primal_step_size);
  }
  else {
    cost += dtau * constraints_->costSlackBarrier(constraints_data_);
  }
  return cost;
}


inline double TerminalUnParNMPC::constraintViolation(
    Robot& robot, const double t, const double dtau, 
    const Eigen::VectorXd& q_prev, const Eigen::VectorXd& v_prev, 
    const SplitSolution& s) {
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  stateequation::ComputeBackwardEulerResidual(robot, dtau, q_prev, v_prev, s, 
                                              kkt_residual_);
  unconstrained_dynamics_.computeUnconstrainedDynamicsResidual(robot, s);
  double violation = 0;
  violation += stateequation::L1NormStateEuqationResidual(kkt_residual_);
  violation += unconstrained_dynamics_.l1NormUnconstrainedDynamicsResidual(dtau);
  violation += dtau * constraints_->l1NormPrimalResidual(constraints_data_);
  return violation;
}


template <typename MatrixType>
inline void TerminalUnParNMPC::computeTerminalCostHessian(
    Robot& robot, const double t, const SplitSolution& s, 
    const Eigen::MatrixBase<MatrixType>& Qxx) {
  assert(Qxx.rows() == 2*robot.dimv());
  assert(Qxx.cols() == 2*robot.dimv());
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  kkt_matrix_.Qxx().setZero();
  cost_->computeTerminalCostHessian(robot, cost_data_, t, s, kkt_matrix_);
  const_cast<Eigen::MatrixBase<MatrixType>&>(Qxx) = kkt_matrix_.Qxx();
}

} // namespace idocp

#endif // IDOCP_TERMINAL_UNPARNMPC_HXX_ 