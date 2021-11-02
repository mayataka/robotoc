#ifndef ROBOTOC_TERMINAL_UNCONSTR_PARNMPC_HXX_
#define ROBOTOC_TERMINAL_UNCONSTR_PARNMPC_HXX_

#include "robotoc/unconstr/terminal_unconstr_parnmpc.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>

namespace robotoc {

inline TerminalUnconstrParNMPC::TerminalUnconstrParNMPC(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot, 0)),
    unconstr_dynamics_(robot),
    use_kinematics_(false),
    stage_cost_(0) {
  if (cost_->useKinematics() || constraints_->useKinematics()) {
    use_kinematics_ = true;
  }
  try {
    if (robot.hasFloatingBase()) {
      throw std::logic_error(
          "robot has floating base: robot should have no constraints!");
    }
    if (robot.maxPointContacts() > 0) {
      throw std::logic_error(
          "robot can have contacts: robot should have no constraints!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline TerminalUnconstrParNMPC::TerminalUnconstrParNMPC() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    unconstr_dynamics_(),
    use_kinematics_(false),
    stage_cost_(0) {
}


inline TerminalUnconstrParNMPC::~TerminalUnconstrParNMPC() {
}


inline bool TerminalUnconstrParNMPC::isFeasible(Robot& robot, 
                                                const SplitSolution& s) {
  return constraints_->isFeasible(robot, constraints_data_, s);
}


inline void TerminalUnconstrParNMPC::initConstraints(Robot& robot, 
                                                     const int time_step, 
                                                     const SplitSolution& s) { 
  assert(time_step >= 0);
  constraints_data_ = constraints_->createConstraintsData(robot, time_step);
  constraints_->setSlackAndDual(robot, constraints_data_, s);
}


inline void TerminalUnconstrParNMPC::evalOCP(Robot& robot, const double t, 
                                             const double dt, 
                                             const Eigen::VectorXd& q_prev, 
                                             const Eigen::VectorXd& v_prev, 
                                             const SplitSolution& s, 
                                             SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  if (use_kinematics_) {
    robot.updateKinematics(s.q);
  }
  kkt_residual.setZero();
  stage_cost_  = cost_->evalStageCost(robot, cost_data_, t, dt, s);
  stage_cost_ += cost_->evalTerminalCost(robot, cost_data_, t, s);
  constraints_->evalConstraint(robot, constraints_data_, s);
  stage_cost_ += constraints_data_.logBarrier();
  unconstr::stateequation::computeBackwardEulerResidual(dt, q_prev, v_prev, s, 
                                                        kkt_residual);
  unconstr_dynamics_.evalUnconstrDynamics(robot, s);
}


inline void TerminalUnconstrParNMPC::computeKKTResidual(
    Robot& robot, const double t, const double dt, 
    const Eigen::VectorXd& q_prev, const Eigen::VectorXd& v_prev, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  if (use_kinematics_) {
    robot.updateKinematics(s.q);
  }
  kkt_residual.setZero();
  stage_cost_  = cost_->linearizeStageCost(robot, cost_data_, t, dt, s, 
                                           kkt_residual);
  stage_cost_ += cost_->linearizeTerminalCost(robot, cost_data_, t, s, 
                                              kkt_residual);
  constraints_->linearizeConstraints(robot, constraints_data_, s, kkt_residual);
  stage_cost_ += constraints_data_.logBarrier();
  unconstr::stateequation::linearizeBackwardEulerTerminal(dt, q_prev, v_prev, s,  
                                                          kkt_matrix, kkt_residual);
  unconstr_dynamics_.linearizeUnconstrDynamics(robot, dt, s, kkt_residual);
  kkt_residual.kkt_error = KKTError(kkt_residual, dt);
}


inline void TerminalUnconstrParNMPC::computeKKTSystem(
    Robot& robot, const double t, const double dt, 
    const Eigen::VectorXd& q_prev, const Eigen::VectorXd& v_prev, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  if (use_kinematics_) {
    robot.updateKinematics(s.q);
  }
  kkt_matrix.setZero();
  kkt_residual.setZero();
  stage_cost_  = cost_->quadratizeStageCost(robot, cost_data_, t, dt, s, 
                                            kkt_residual, kkt_matrix);
  stage_cost_ += cost_->quadratizeTerminalCost(robot, cost_data_, t, s, 
                                               kkt_residual, kkt_matrix);
  constraints_->linearizeConstraints(robot, constraints_data_, s, kkt_residual);
  stage_cost_ += constraints_data_.logBarrier();
  unconstr::stateequation::linearizeBackwardEulerTerminal(dt, q_prev, v_prev, s,  
                                                          kkt_matrix, kkt_residual);
  unconstr_dynamics_.linearizeUnconstrDynamics(robot, dt, s, kkt_residual);
  kkt_residual.kkt_error = KKTError(kkt_residual, dt);
  constraints_->condenseSlackAndDual(constraints_data_, s, kkt_matrix, 
                                     kkt_residual);
  unconstr_dynamics_.condenseUnconstrDynamics(kkt_matrix, kkt_residual);
}


inline void TerminalUnconstrParNMPC::expandPrimalAndDual(
    const double dt, const SplitSolution& s, const SplitKKTMatrix& kkt_matrix, 
    const SplitKKTResidual& kkt_residual, SplitDirection& d) {
  assert(dt > 0);
  unconstr_dynamics_.expandPrimal(d);
  unconstr_dynamics_.expandDual(dt, kkt_matrix, kkt_residual, d);
  constraints_->expandSlackAndDual(constraints_data_, s, d);
}


inline double TerminalUnconstrParNMPC::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


inline double TerminalUnconstrParNMPC::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


inline void TerminalUnconstrParNMPC::updatePrimal(const Robot& robot, 
                                                  const double primal_step_size, 
                                                  const SplitDirection& d, 
                                                  SplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  constraints_->updateSlack(constraints_data_, primal_step_size);
}


inline void TerminalUnconstrParNMPC::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(constraints_data_, dual_step_size);
}


inline double TerminalUnconstrParNMPC::KKTError(
    const SplitKKTResidual& kkt_residual, const double dt) const {
  assert(dt > 0);
  double err = 0;
  err += kkt_residual.KKTError();
  err += (dt*dt) * unconstr_dynamics_.KKTError();
  err += constraints_data_.KKTError();
  return err;
}


inline double TerminalUnconstrParNMPC::stageCost() const {
  return stage_cost_;
}


inline double TerminalUnconstrParNMPC::constraintViolation(
    const SplitKKTResidual& kkt_residual, const double dt) const {
  double vio = 0;
  vio += kkt_residual.constraintViolation();
  vio += dt * unconstr_dynamics_.constraintViolation();
  vio += constraints_data_.constraintViolation();
  return vio;
}


inline void TerminalUnconstrParNMPC::evalTerminalCostHessian(
    Robot& robot, const double t, const SplitSolution& s, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) {
  if (use_kinematics_) {
    robot.updateKinematics(s.q);
  }
  kkt_matrix.setZero();
  kkt_residual.setZero();
  double cost = cost_->quadratizeTerminalCost(robot, cost_data_, t, s, 
                                              kkt_residual, kkt_matrix);
}

} // namespace robotoc

#endif // ROBOTOC_TERMINAL_UNCONSTR_PARNMPC_HXX_ 