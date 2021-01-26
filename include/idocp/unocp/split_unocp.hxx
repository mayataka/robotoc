#ifndef IDOCP_SPLIT_UNOCP_HXX_
#define IDOCP_SPLIT_UNOCP_HXX_

#include "idocp/unocp/split_unocp.hpp"

#include <stdexcept>
#include <cassert>

namespace idocp {

inline SplitUnOCP::SplitUnOCP(const Robot& robot, 
                              const std::shared_ptr<CostFunction>& cost,
                              const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot, 0)),
    unconstrained_dynamics_(robot),
    use_kinematics_(false),
    kkt_matrix_(robot),
    kkt_residual_(robot) {
  if (cost_->useKinematics() || constraints_->useKinematics()) {
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


inline SplitUnOCP::SplitUnOCP() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    unconstrained_dynamics_(),
    use_kinematics_(false),
    kkt_matrix_(),
    kkt_residual_() {
}


inline SplitUnOCP::~SplitUnOCP() {
}


inline bool SplitUnOCP::isFeasible(Robot& robot, const SplitSolution& s) {
  return constraints_->isFeasible(robot, constraints_data_, s);
}


inline void SplitUnOCP::initConstraints(Robot& robot, const int time_step, 
                                        const SplitSolution& s) { 
  assert(time_step >= 0);
  constraints_data_ = constraints_->createConstraintsData(robot, time_step);
  constraints_->setSlackAndDual(robot, constraints_data_, s);
}


inline void SplitUnOCP::linearizeOCP(Robot& robot, const double t, 
                                     const double dtau, 
                                     const Eigen::VectorXd& q_prev, 
                                     const SplitSolution& s, 
                                     const SplitSolution& s_next, 
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
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s,
                                    kkt_residual_);
  stateequation::linearizeForwardEuler(robot, dtau, q_prev, s, s_next, 
                                       kkt_matrix_, kkt_residual_);
  unconstrained_dynamics_.linearizeUnconstrainedDynamics(robot, dtau, s, 
                                                         kkt_residual_);
  cost_->computeStageCostHessian(robot, cost_data_, t, dtau, s, kkt_matrix_);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dtau, s, 
                                     kkt_matrix_, kkt_residual_);
  unconstrained_dynamics_.condenseUnconstrainedDynamics(
      kkt_matrix_, kkt_residual_, unkkt_matrix, unkkt_residual);
}


inline void SplitUnOCP::computeCondensedDirection(Robot& robot, 
                                                  const double dtau, 
                                                  const SplitSolution& s, 
                                                  SplitDirection& d) {
  unconstrained_dynamics_.computeCondensedDirection(dtau, kkt_matrix_, 
                                                    kkt_residual_, d);
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, s, d);
}


inline double SplitUnOCP::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


inline double SplitUnOCP::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


inline void SplitUnOCP::updatePrimal(const Robot& robot, 
                                     const double primal_step_size, 
                                     const SplitDirection& d, 
                                     SplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  constraints_->updateSlack(constraints_data_, primal_step_size);
}


inline void SplitUnOCP::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(constraints_data_, dual_step_size);
}


inline void SplitUnOCP::computeKKTResidual(Robot& robot, const double t, 
                                           const double dtau, 
                                           const Eigen::VectorXd& q_prev, 
                                           const SplitSolution& s,
                                           const SplitSolution& s_next) {
  assert(dtau >= 0);
  assert(q_prev.size() == robot.dimq());
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  kkt_residual_.setZero();
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s,
                                    kkt_residual_);
  stateequation::linearizeForwardEuler(robot, dtau, q_prev, s, s_next, 
                                       kkt_matrix_, kkt_residual_);
  unconstrained_dynamics_.linearizeUnconstrainedDynamics(robot, dtau, s, 
                                                         kkt_residual_);
}


inline double SplitUnOCP::squaredNormKKTResidual(const double dtau) const {
  double error = 0;
  error += kkt_residual_.lx().squaredNorm();
  error += kkt_residual_.la.squaredNorm();
  error += kkt_residual_.lu().squaredNorm();
  error += stateequation::squaredNormStateEuqationResidual(kkt_residual_);
  error += unconstrained_dynamics_.squaredNormUnconstrainedDynamicsResidual(dtau);
  error += dtau * dtau * constraints_->squaredNormPrimalAndDualResidual(constraints_data_);
  return error;
}


inline double SplitUnOCP::stageCost(Robot& robot, const double t,  
                                    const double dtau, const SplitSolution& s, 
                                    const double primal_step_size) {
  assert(dtau >= 0);
  assert(primal_step_size >= 0);
  assert(primal_step_size <= 1);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  double cost = 0;
  cost += cost_->computeStageCost(robot, cost_data_, t, dtau, s);
  if (primal_step_size > 0) {
    cost += dtau * constraints_->costSlackBarrier(constraints_data_, 
                                                  primal_step_size);
  }
  else {
    cost += dtau * constraints_->costSlackBarrier(constraints_data_);
  }
  return cost;
}


inline double SplitUnOCP::constraintViolation(Robot& robot, const double t, 
                                              const double dtau, 
                                              const SplitSolution& s, 
                                              const Eigen::VectorXd& q_next, 
                                              const Eigen::VectorXd& v_next) {
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  stateequation::computeForwardEulerResidual(robot, dtau, s, q_next, v_next, 
                                             kkt_residual_);
  unconstrained_dynamics_.computeUnconstrainedDynamicsResidual(robot, s);
  double violation = 0;
  violation += stateequation::l1NormStateEuqationResidual(kkt_residual_);
  violation += unconstrained_dynamics_.l1NormUnconstrainedDynamicsResidual(dtau);
  violation += dtau * constraints_->l1NormPrimalResidual(constraints_data_);
  return violation;
}

} // namespace idocp

#endif // IDOCP_SPLIT_UNOCP_HXX_