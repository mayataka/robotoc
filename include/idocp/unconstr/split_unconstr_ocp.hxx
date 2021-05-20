#ifndef IDOCP_SPLIT_UNCONSTR_OCP_HXX_
#define IDOCP_SPLIT_UNCONSTR_OCP_HXX_

#include "idocp/unocp/split_unconstr_ocp.hpp"

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
    kkt_matrix_(robot),
    kkt_residual_(robot) {
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
    kkt_matrix_(),
    kkt_residual_() {
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
                                           SplitUnKKTMatrix& unkkt_matrix,
                                           SplitUnKKTResidual& unkkt_residual) {
  assert(dt > 0);
  assert(q_prev.size() == robot.dimq());
  if (use_kinematics_) {
    robot.updateKinematics(s.q);
  }
  kkt_matrix_.setZero();
  kkt_residual_.setZero();
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dt, s, kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, dt, s, kkt_residual_);
  stateequation::linearizeForwardEuler(robot, dt, q_prev, s, s_next, 
                                       kkt_matrix_, kkt_residual_);
  unconstr_dynamics_.linearizeUnconstrDynamics(robot, dt, s, kkt_residual_);
  cost_->computeStageCostHessian(robot, cost_data_, t, dt, s, kkt_matrix_);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dt, s, 
                                     kkt_matrix_, kkt_residual_);
  unconstr_dynamics_.condenseUnconstrDynamics(kkt_matrix_, kkt_residual_, 
                                              unkkt_matrix, unkkt_residual);
}


inline void SplitUnconstrOCP::computeCondensedDirection(Robot& robot, 
                                                        const double dt, 
                                                        const SplitSolution& s, 
                                                        SplitDirection& d) {
  assert(dt > 0);
  unconstr_dynamics_.computeCondensedDirection(dt, kkt_matrix_, kkt_residual_, d);
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
                                                 const SplitSolution& s_next) {
  assert(dt > 0);
  assert(q_prev.size() == robot.dimq());
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  kkt_residual_.setZero();
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dt, s, kkt_residual_);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  constraints_->augmentDualResidual(robot, constraints_data_, dt, s, kkt_residual_);
  stateequation::linearizeForwardEuler(robot, dt, q_prev, s, s_next, 
                                       kkt_matrix_, kkt_residual_);
  unconstr_dynamics_.linearizeUnconstDynamics(robot, dt, s, kkt_residual_);
}


inline double SplitUnconstrOCP::squaredNormKKTResidual(const double dt) const {
  assert(dt > 0);
  double error = 0;
  error += kkt_residual_.lx.squaredNorm();
  error += kkt_residual_.la.squaredNorm();
  error += kkt_residual_.lu.squaredNorm();
  error += stateequation::squaredNormStateEuqationResidual(kkt_residual_);
  error += unconstr_dynamics_.squaredNormUnconstrDynamicsResidual(dt);
  error += dt * dt * constraints_->squaredNormPrimalAndDualResidual(constraints_data_);
  return error;
}


inline double SplitUnconstrOCP::stageCost(Robot& robot, const double t,  
                                    const double dt, const SplitSolution& s, 
                                    const double primal_step_size) {
  assert(dt > 0);
  assert(primal_step_size >= 0);
  assert(primal_step_size <= 1);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
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


inline double SplitUnconstrOCP::constraintViolation(Robot& robot, const double t, 
                                              const double dt, 
                                              const SplitSolution& s, 
                                              const Eigen::VectorXd& q_next, 
                                              const Eigen::VectorXd& v_next) {
  assert(dt > 0);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  stateequation::computeForwardEulerResidual(robot, dt, s, q_next, v_next, 
                                             kkt_residual_);
  unconstr_dynamics_.computeUnconstrDynamicsResidual(robot, s);
  double violation = 0;
  violation += stateequation::l1NormStateEuqationResidual(kkt_residual_);
  violation += unconstr_dynamics_.l1NormUnconstrDynamicsResidual(dt);
  violation += dt * constraints_->l1NormPrimalResidual(constraints_data_);
  return violation;
}

} // namespace idocp

#endif // IDOCP_SPLIT_UNCONSTR_OCP_HXX_ 