#ifndef IDOCP_SPLIT_IMPULSE_OCP_HXX_
#define IDOCP_SPLIT_IMPULSE_OCP_HXX_

#include "idocp/impulse/split_impulse_ocp.hpp"

#include <cassert>

namespace idocp {

inline SplitImpulseOCP::SplitImpulseOCP(
    const Robot& robot, const std::shared_ptr<ImpulseCostFunction>& cost, 
    const std::shared_ptr<ImpulseConstraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot)),
    impulse_dynamics_(robot) {
}


inline SplitImpulseOCP::SplitImpulseOCP() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    impulse_dynamics_() {
}


inline SplitImpulseOCP::~SplitImpulseOCP() {
}


inline bool SplitImpulseOCP::isFeasible(Robot& robot, 
                                        const ImpulseSplitSolution& s) {
  return constraints_->isFeasible(robot, constraints_data_, s);
}


inline void SplitImpulseOCP::initConstraints(Robot& robot,
                                             const ImpulseSplitSolution& s) { 
  constraints_data_ = constraints_->createConstraintsData(robot);
  constraints_->setSlackAndDual(robot, constraints_data_, s);
}


inline void SplitImpulseOCP::linearizeOCP(
    Robot& robot, const ImpulseStatus& impulse_status, const double t,  
    const Eigen::VectorXd& q_prev, const ImpulseSplitSolution& s, 
    const SplitSolution& s_next, ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual, 
    const bool is_state_constraint_valid) {
  kkt_matrix.setImpulseStatus(impulse_status);
  kkt_residual.setImpulseStatus(impulse_status);
  robot.updateKinematics(s.q, s.v+s.dv);
  // condensing the impulse dynamics
  kkt_residual.setZero();
  kkt_matrix.setZero();
  cost_->computeStageCostDerivatives(robot, cost_data_, t, s, kkt_residual);
  constraints_->augmentDualResidual(robot, constraints_data_, s, kkt_residual);
  stateequation::LinearizeImpulseForwardEuler(robot, q_prev, s, s_next, 
                                              kkt_matrix, kkt_residual);
  impulse_dynamics_.linearizeImpulseDynamics(robot, impulse_status, s,
                                             kkt_matrix, kkt_residual,
                                             is_state_constraint_valid);
  cost_->computeStageCostHessian(robot, cost_data_, t, s, kkt_matrix);
  constraints_->condenseSlackAndDual(robot, constraints_data_, s, kkt_matrix, 
                                     kkt_residual);
  impulse_dynamics_.condenseImpulseDynamics(robot, impulse_status, 
                                            kkt_matrix, kkt_residual);
}


inline void SplitImpulseOCP::computeCondensedPrimalDirection(
    Robot& robot, const ImpulseSplitSolution& s, ImpulseSplitDirection& d) {
  impulse_dynamics_.computeCondensedPrimalDirection(robot, d);
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, s, d);
}


inline void SplitImpulseOCP::computeCondensedDualDirection(
    const Robot& robot, const ImpulseSplitKKTMatrix& kkt_matrix, 
    const ImpulseSplitKKTResidual& kkt_residual, const SplitDirection& d_next, 
    ImpulseSplitDirection& d) {
  impulse_dynamics_.computeCondensedDualDirection(robot, kkt_matrix,
                                                  kkt_residual, 
                                                  d_next.dgmm(), d);
}


inline double SplitImpulseOCP::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


inline double SplitImpulseOCP::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


inline void SplitImpulseOCP::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(constraints_data_, dual_step_size);
}


inline void SplitImpulseOCP::updatePrimal(
    const Robot& robot, const double primal_step_size, 
    const ImpulseSplitDirection& d, ImpulseSplitSolution& s, 
    const bool is_state_constraint_valid) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d, is_state_constraint_valid);
  constraints_->updateSlack(constraints_data_, primal_step_size);
}


inline void SplitImpulseOCP::computeKKTResidual(
    Robot& robot, const ImpulseStatus& impulse_status, const double t, 
    const Eigen::VectorXd& q_prev, const ImpulseSplitSolution& s, 
    const SplitSolution& s_next, ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual, 
    const bool is_state_constraint_valid) {
  assert(q_prev.size() == robot.dimq());
  kkt_matrix.setImpulseStatus(impulse_status);
  kkt_residual.setImpulseStatus(impulse_status);
  kkt_residual.setZero();
  robot.updateKinematics(s.q, s.v+s.dv);
  cost_->computeStageCostDerivatives(robot, cost_data_, t, s, kkt_residual);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  constraints_->augmentDualResidual(robot, constraints_data_, s, kkt_residual);
  stateequation::LinearizeImpulseForwardEuler(robot, q_prev, s, s_next, 
                                              kkt_matrix, kkt_residual);
  impulse_dynamics_.linearizeImpulseDynamics(robot, impulse_status, s, 
                                             kkt_matrix, kkt_residual, 
                                             is_state_constraint_valid);
}


inline double SplitImpulseOCP::squaredNormKKTResidual(
    const ImpulseSplitKKTResidual& kkt_residual, 
    const bool is_state_constraint_valid) const {
  double error = 0;
  error += kkt_residual.lx().squaredNorm();
  error += kkt_residual.ldv.squaredNorm();
  error += kkt_residual.lf().squaredNorm();
  error += stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  error += impulse_dynamics_.squaredNormImpulseDynamicsResidual(
      kkt_residual, is_state_constraint_valid);
  error += constraints_->squaredNormPrimalAndDualResidual(constraints_data_);
  return error;
}


inline double SplitImpulseOCP::stageCost(Robot& robot, const double t, 
                                         const ImpulseSplitSolution& s, 
                                         const double primal_step_size) {
  assert(primal_step_size >= 0);
  assert(primal_step_size <= 1);
  robot.updateKinematics(s.q, s.v+s.dv);
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, s);
  if (primal_step_size > 0) {
    cost += constraints_->costSlackBarrier(constraints_data_, primal_step_size);
  }
  else {
    cost += constraints_->costSlackBarrier(constraints_data_);
  }
  return cost;
}


inline double SplitImpulseOCP::constraintViolation(
    Robot& robot, const ImpulseStatus& impulse_status, const double t, 
    const ImpulseSplitSolution& s, const Eigen::VectorXd& q_next, 
    const Eigen::VectorXd& v_next, ImpulseSplitKKTResidual& kkt_residual,
    const bool is_state_constraint_valid) {
  kkt_residual.setImpulseStatus(impulse_status);
  kkt_residual.setZero();
  robot.updateKinematics(s.q, s.v+s.dv);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  stateequation::ComputeImpulseForwardEulerResidual(robot, s, q_next, v_next, 
                                                    kkt_residual);
  impulse_dynamics_.computeImpulseDynamicsResidual(robot, impulse_status, s, 
                                                   kkt_residual,
                                                   is_state_constraint_valid);
  double violation = 0;
  violation += constraints_->l1NormPrimalResidual(constraints_data_);
  violation += stateequation::L1NormStateEuqationResidual(kkt_residual);
  violation += impulse_dynamics_.l1NormImpulseDynamicsResidual(
      kkt_residual, is_state_constraint_valid);
  return violation;
}

} // namespace idocp

#endif // IDOCP_SPLIT_IMPULSE_OCP_HXX_ 