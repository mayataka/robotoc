#ifndef ROBOTOC_SPLIT_UNCONSTR_PARNMPC_HXX_
#define ROBOTOC_SPLIT_UNCONSTR_PARNMPC_HXX_

#include "robotoc/unconstr/split_unconstr_parnmpc.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>

namespace robotoc {

inline SplitUnconstrParNMPC::SplitUnconstrParNMPC(
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


inline SplitUnconstrParNMPC::SplitUnconstrParNMPC() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    unconstr_dynamics_(),
    use_kinematics_(false),
    stage_cost_(0) {
}


inline SplitUnconstrParNMPC::~SplitUnconstrParNMPC() {
}


inline bool SplitUnconstrParNMPC::isFeasible(Robot& robot, 
                                             const SplitSolution& s) {
  return constraints_->isFeasible(robot, constraints_data_, s);
}


inline void SplitUnconstrParNMPC::initConstraints(Robot& robot, 
                                                  const int time_step, 
                                                  const SplitSolution& s) { 
  assert(time_step >= 0);
  constraints_data_ = constraints_->createConstraintsData(robot, time_step);
  constraints_->setSlackAndDual(robot, constraints_data_, s);
}


inline void SplitUnconstrParNMPC::evalOCP(Robot& robot, const double t, 
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
  stage_cost_ = cost_->computeStageCost(robot, cost_data_, t, dt, s);
  constraints_->evalConstraint(robot, constraints_data_, s);
  stage_cost_ += dt * constraints_data_.logBarrier();
  unconstr::stateequation::computeBackwardEulerResidual(dt, q_prev, v_prev, s, 
                                                        kkt_residual);
  unconstr_dynamics_.computeUnconstrDynamicsResidual(robot, s);
}


inline void SplitUnconstrParNMPC::computeKKTResidual(Robot& robot, const double t, 
                                                     const double dt, 
                                                     const Eigen::VectorXd& q_prev, 
                                                     const Eigen::VectorXd& v_prev, 
                                                     const SplitSolution& s, 
                                                     const SplitSolution& s_next, 
                                                     SplitKKTMatrix& kkt_matrix, 
                                                     SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  if (use_kinematics_) {
    robot.updateKinematics(s.q);
  }
  kkt_residual.setZero();
  stage_cost_ = cost_->linearizeStageCost(robot, cost_data_, t, dt, s, 
                                          kkt_residual);
  constraints_->linearizeConstraints(robot, constraints_data_, dt, s, 
                                     kkt_residual);
  stage_cost_ += dt * constraints_data_.logBarrier();
  unconstr::stateequation::linearizeBackwardEuler(dt, q_prev, v_prev, s, s_next, 
                                                  kkt_matrix, kkt_residual);
  unconstr_dynamics_.linearizeUnconstrDynamics(robot, dt, s, kkt_residual);
}


inline void SplitUnconstrParNMPC::computeKKTSystem(Robot& robot, const double t, 
                                                   const double dt, 
                                                   const Eigen::VectorXd& q_prev, 
                                                   const Eigen::VectorXd& v_prev, 
                                                   const SplitSolution& s, 
                                                   const SplitSolution& s_next, 
                                                   SplitKKTMatrix& kkt_matrix,
                                                   SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  if (use_kinematics_) {
    robot.updateKinematics(s.q);
  }
  kkt_matrix.setZero();
  kkt_residual.setZero();
  stage_cost_ = cost_->quadratizeStageCost(robot, cost_data_, t, dt, s, 
                                           kkt_residual, kkt_matrix);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dt, s, 
                                     kkt_matrix, kkt_residual);
  stage_cost_ += dt * constraints_data_.logBarrier();
  unconstr::stateequation::linearizeBackwardEuler(dt, q_prev, v_prev, s, s_next, 
                                                  kkt_matrix, kkt_residual);
  unconstr_dynamics_.linearizeUnconstrDynamics(robot, dt, s, kkt_residual);
  unconstr_dynamics_.condenseUnconstrDynamics(kkt_matrix, kkt_residual);
}


inline void SplitUnconstrParNMPC::expandPrimalAndDual(
    const double dt, const SplitSolution& s, const SplitKKTMatrix& kkt_matrix, 
    const SplitKKTResidual& kkt_residual, SplitDirection& d) {
  assert(dt > 0);
  unconstr_dynamics_.expandPrimal(d);
  unconstr_dynamics_.expandDual(dt, kkt_matrix, kkt_residual, d);
  constraints_->expandSlackAndDual(constraints_data_, s, d);
}


inline double SplitUnconstrParNMPC::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


inline double SplitUnconstrParNMPC::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


inline void SplitUnconstrParNMPC::updatePrimal(const Robot& robot, 
                                               const double primal_step_size, 
                                               const SplitDirection& d, 
                                               SplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  constraints_->updateSlack(constraints_data_, primal_step_size);
}


inline void SplitUnconstrParNMPC::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(constraints_data_, dual_step_size);
}


inline double SplitUnconstrParNMPC::KKTError(
    const SplitKKTResidual& kkt_residual, const double dt) const {
  assert(dt > 0);
  double err = 0;
  err += kkt_residual.KKTError();
  err += (dt*dt) * unconstr_dynamics_.KKTError();
  err += (dt*dt) * constraints_data_.KKTError();
  return err;
}


inline double SplitUnconstrParNMPC::stageCost() const {
  return stage_cost_;
}


inline double SplitUnconstrParNMPC::constraintViolation(
    const SplitKKTResidual& kkt_residual, const double dt) const {
  double vio = 0;
  vio += kkt_residual.constraintViolation();
  vio += dt * unconstr_dynamics_.constraintViolation();
  vio += dt * constraints_data_.constraintViolation();
  return vio;
}

} // namespace robotoc

#endif // ROBOTOC_SPLIT_UNCONSTR_PARNMPC_HXX_ 