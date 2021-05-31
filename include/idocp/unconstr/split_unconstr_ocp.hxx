#ifndef IDOCP_SPLIT_UNCONSTR_OCP_HXX_
#define IDOCP_SPLIT_UNCONSTR_OCP_HXX_

#include "idocp/unconstr/split_unconstr_ocp.hpp"

#include <stdexcept>
#include <cassert>

namespace idocp {

inline SplitUnconstrOCP::SplitUnconstrOCP(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost,
    const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot, 0)),
    unconstr_dynamics_(robot),
    use_kinematics_(false),
    stage_cost_(0),
    constraint_violation_(0) {
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


inline SplitUnconstrOCP::SplitUnconstrOCP() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    unconstr_dynamics_(),
    use_kinematics_(false),
    stage_cost_(0),
    constraint_violation_(0) {
}


inline SplitUnconstrOCP::~SplitUnconstrOCP() {
}


inline bool SplitUnconstrOCP::isFeasible(Robot& robot, const SplitSolution& s) {
  return constraints_->isFeasible(robot, constraints_data_, s);
}


inline void SplitUnconstrOCP::initConstraints(Robot& robot, const int time_step, 
                                              const SplitSolution& s) { 
  assert(time_step >= 0);
  constraints_data_ = constraints_->createConstraintsData(robot, time_step);
  constraints_->setSlackAndDual(robot, constraints_data_, s);
}


inline void SplitUnconstrOCP::linearizeOCP(Robot& robot, const double t, 
                                           const double dt, 
                                           const Eigen::VectorXd& q_prev, 
                                           const SplitSolution& s, 
                                           const SplitSolution& s_next, 
                                           SplitKKTMatrix& kkt_matrix,
                                           SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  assert(q_prev.size() == robot.dimq());
  if (use_kinematics_) {
    robot.updateKinematics(s.q);
  }
  kkt_matrix.setZero();
  kkt_residual.setZero();
  stage_cost_ = cost_->quadratizeStageCost(robot, cost_data_, t, dt, s, 
                                           kkt_residual, kkt_matrix);
  constraints_->augmentDualResidual(robot, constraints_data_, dt, s, kkt_residual);
  stateequation::linearizeForwardEuler(robot, dt, q_prev, s, s_next, 
                                       kkt_matrix, kkt_residual);
  unconstr_dynamics_.linearizeUnconstrDynamics(robot, dt, s, kkt_residual);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dt, s, 
                                     kkt_matrix, kkt_residual);
  unconstr_dynamics_.condenseUnconstrDynamics(kkt_matrix, kkt_residual);
}


inline void SplitUnconstrOCP::computeCondensedDirection(
    Robot& robot, const double dt, const SplitSolution& s, 
    const SplitKKTMatrix& kkt_matrix, const SplitKKTResidual& kkt_residual,
    SplitDirection& d) {
  assert(dt > 0);
  unconstr_dynamics_.computeCondensedDirection(dt, kkt_matrix, kkt_residual, d);
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, s, d);
}


inline double SplitUnconstrOCP::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


inline double SplitUnconstrOCP::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


inline void SplitUnconstrOCP::updatePrimal(const Robot& robot, 
                                           const double primal_step_size, 
                                           const SplitDirection& d, 
                                           SplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  constraints_->updateSlack(constraints_data_, primal_step_size);
}


inline void SplitUnconstrOCP::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(constraints_data_, dual_step_size);
}


inline void SplitUnconstrOCP::computeKKTResidual(Robot& robot, const double t, 
                                                 const double dt, 
                                                 const Eigen::VectorXd& q_prev, 
                                                 const SplitSolution& s,
                                                 const SplitSolution& s_next,
                                                 SplitKKTMatrix& kkt_matrix, 
                                                 SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  assert(q_prev.size() == robot.dimq());
  if (use_kinematics_) {
    robot.updateKinematics(s.q);
  }
  kkt_residual.setZero();
  stage_cost_ = cost_->linearizeStageCost(robot, cost_data_, t, dt, s, 
                                          kkt_residual);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  constraints_->augmentDualResidual(robot, constraints_data_, dt, s, kkt_residual);
  stateequation::linearizeForwardEuler(robot, dt, q_prev, s, s_next, 
                                       kkt_matrix, kkt_residual);
  unconstr_dynamics_.linearizeUnconstrDynamics(robot, dt, s, kkt_residual);
}


inline double SplitUnconstrOCP::squaredNormKKTResidual(
    const SplitKKTResidual& kkt_residual, const double dt) const {
  assert(dt > 0);
  double error = 0;
  error += kkt_residual.lx.squaredNorm();
  error += kkt_residual.la.squaredNorm();
  error += kkt_residual.lu.squaredNorm();
  error += stateequation::squaredNormStateEuqationResidual(kkt_residual);
  error += unconstr_dynamics_.squaredNormUnconstrDynamicsResidual(dt);
  error += dt * dt * constraints_->squaredNormPrimalAndDualResidual(constraints_data_);
  return error;
}


inline double SplitUnconstrOCP::stageCost(Robot& robot, const double t,  
                                          const double dt, 
                                          const SplitSolution& s, 
                                          const double primal_step_size) {
  assert(dt > 0);
  assert(primal_step_size >= 0);
  assert(primal_step_size <= 1);
  if (use_kinematics_) {
    robot.updateKinematics(s.q);
  }
  double cost = 0;
  cost += cost_->computeStageCost(robot, cost_data_, t, dt, s);
  if (primal_step_size > 0) {
    cost += dt * constraints_->costSlackBarrier(constraints_data_, 
                                                primal_step_size);
  }
  else {
    cost += dt * constraints_->costSlackBarrier(constraints_data_);
  }
  return cost;
}


inline double SplitUnconstrOCP::constraintViolation(
    Robot& robot, const double t, const double dt, const SplitSolution& s, 
    const Eigen::VectorXd& q_next, const Eigen::VectorXd& v_next,
    SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  if (use_kinematics_) {
    robot.updateKinematics(s.q);
  }
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  stateequation::computeForwardEulerResidual(robot, dt, s, q_next, v_next, 
                                             kkt_residual);
  unconstr_dynamics_.computeUnconstrDynamicsResidual(robot, s);
  double violation = 0;
  violation += stateequation::l1NormStateEuqationResidual(kkt_residual);
  violation += unconstr_dynamics_.l1NormUnconstrDynamicsResidual(dt);
  violation += dt * constraints_->l1NormPrimalResidual(constraints_data_);
  return violation;
}

} // namespace idocp

#endif // IDOCP_SPLIT_UNCONSTR_OCP_HXX_ 