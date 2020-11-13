#include "idocp/impulse/split_impulse_ocp.hpp"

#include <cassert>


namespace idocp {

SplitImpulseOCP::SplitImpulseOCP(
    const Robot& robot, const std::shared_ptr<ImpulseCostFunction>& cost, 
    const std::shared_ptr<ImpulseConstraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot)),
    kkt_residual_(robot),
    kkt_matrix_(robot),
    impulse_dynamics_(robot),
    riccati_factorizer_(robot) {
}


SplitImpulseOCP::SplitImpulseOCP() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    kkt_residual_(),
    kkt_matrix_(),
    impulse_dynamics_(),
    riccati_factorizer_() {
}


SplitImpulseOCP::~SplitImpulseOCP() {
}


bool SplitImpulseOCP::isFeasible(Robot& robot, const ImpulseSplitSolution& s) {
  return constraints_->isFeasible(robot, constraints_data_, s);
}


void SplitImpulseOCP::initConstraints(Robot& robot,
                                      const ImpulseSplitSolution& s) { 
  constraints_data_ = constraints_->createConstraintsData(robot);
  constraints_->setSlackAndDual(robot, constraints_data_, s);
}


void SplitImpulseOCP::linearizeOCP(Robot& robot, 
                                   const ImpulseStatus& impulse_status,  
                                   const double t,  
                                   const Eigen::VectorXd& q_prev, 
                                   const ImpulseSplitSolution& s, 
                                   const SplitSolution& s_next) {
  setImpulseStatusForKKT(impulse_status);
  robot.updateKinematics(s.q, s.v);
  // condensing the impulse dynamics
  kkt_residual_.setZero();
  kkt_matrix_.setZero();
  cost_->computeStageCostDerivatives(robot, cost_data_, t, s, kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, s, kkt_residual_);
  stateequation::LinearizeImpulseForwardEuler(robot, q_prev, s, s_next, 
                                              kkt_matrix_, kkt_residual_);
  impulse_dynamics_.linearizeImpulseDynamics(robot, impulse_status, s,
                                             kkt_matrix_, kkt_residual_);
  cost_->computeStageCostHessian(robot, cost_data_, t, s, kkt_matrix_);
  constraints_->condenseSlackAndDual(robot, constraints_data_, s, kkt_matrix_, 
                                     kkt_residual_);
  impulse_dynamics_.condenseImpulseDynamics(robot, impulse_status, 
                                            kkt_matrix_, kkt_residual_);
}

void SplitImpulseOCP::getStateConstraint(
    StateConstraintRiccatiFactorization& state_constraint_factorization) const {
  state_constraint_factorization.E() = 
  state_constraint_factorization.e() = 
  state_constraint_factorization.T_impulse().head(dimv_) = 
}

void SplitImpulseOCP::backwardRiccatiRecursion(
    const RiccatiSolution& riccati_next, RiccatiSolution& riccati) {
  riccati_factorizer_.backwardRiccatiRecursion(riccati_next, kkt_matrix_, 
                                               kkt_residual_, riccati);
}


void SplitImpulseOCP::forwardRiccatiRecursion(ImpulseSplitDirection& d,   
                                              SplitDirection& d_next) {
  riccati_factorizer_.forwardRiccatiRecursion(kkt_matrix_, kkt_residual_, d,
                                              d_next);
}


void SplitImpulseOCP::computeCondensedPrimalDirection(
    Robot& robot, const RiccatiSolution& riccati, const ImpulseSplitSolution& s, 
    ImpulseSplitDirection& d) {
  riccati_factorizer_.computeCostateDirection(riccati, d);
  impulse_dynamics_.computeCondensedPrimalDirection(robot, d);
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, s, d);
}


void SplitImpulseOCP::computeCondensedDualDirection(
    Robot& robot, const SplitDirection& d_next, ImpulseSplitDirection& d) {
  impulse_dynamics_.computeCondensedDualDirection(robot, kkt_matrix_,
                                                  kkt_residual_, 
                                                  d_next.dgmm(), d);
}

 
double SplitImpulseOCP::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


double SplitImpulseOCP::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


void SplitImpulseOCP::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(constraints_data_, dual_step_size);
}


void SplitImpulseOCP::updatePrimal(Robot& robot, const double primal_step_size, 
                                   const ImpulseSplitDirection& d, 
                                   ImpulseSplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  constraints_->updateSlack(constraints_data_, primal_step_size);
}


void SplitImpulseOCP::computeKKTResidual(Robot& robot, 
                                         const ImpulseStatus& impulse_status, 
                                         const double t, 
                                         const Eigen::VectorXd& q_prev, 
                                         const ImpulseSplitSolution& s,
                                         const SplitSolution& s_next) {
  assert(q_prev.size() == robot.dimq());
  setImpulseStatusForKKT(impulse_status);
  kkt_residual_.setZero();
  robot.updateKinematics(s.q, s.v);
  cost_->computeStageCostDerivatives(robot, cost_data_, t, s, kkt_residual_);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  constraints_->augmentDualResidual(robot, constraints_data_, s, kkt_residual_);
  stateequation::LinearizeImpulseForwardEuler(robot, q_prev, s, s_next, 
                                              kkt_matrix_, kkt_residual_);
  impulse_dynamics_.linearizeImpulseDynamics(robot, impulse_status, s, 
                                             kkt_matrix_, kkt_residual_);
}


double SplitImpulseOCP::squaredNormKKTResidual() const {
  double error = 0;
  error += kkt_residual_.lx().squaredNorm();
  error += kkt_residual_.ldv.squaredNorm();
  error += kkt_residual_.lf().squaredNorm();
  error += stateequation::SquaredNormStateEuqationResidual(kkt_residual_);
  error += impulse_dynamics_.squaredNormImpulseDynamicsResidual(kkt_residual_);
  error += constraints_->squaredNormPrimalAndDualResidual(constraints_data_);
  return error;
}


double SplitImpulseOCP::stageCost(Robot& robot, const double t, 
                                  const ImpulseSplitSolution& s, 
                                  const double primal_step_size) {
  assert(primal_step_size >= 0);
  assert(primal_step_size <= 1);
  robot.updateKinematics(s.q, s.v);
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


double SplitImpulseOCP::constraintViolation(Robot& robot, 
                                            const ImpulseStatus& impulse_status, 
                                            const double t, 
                                            const ImpulseSplitSolution& s, 
                                            const Eigen::VectorXd& q_next, 
                                            const Eigen::VectorXd& v_next) {
  setImpulseStatusForKKT(impulse_status);
  robot.updateKinematics(s.q, s.v);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  stateequation::ComputeImpulseForwardEulerResidual(robot, s, q_next, v_next, 
                                                    kkt_residual_);
  impulse_dynamics_.computeImpulseDynamicsResidual(robot, impulse_status, s, 
                                                   kkt_residual_);
  double violation = 0;
  violation += constraints_->l1NormPrimalResidual(constraints_data_);
  violation += stateequation::L1NormStateEuqationResidual(kkt_residual_);
  violation += impulse_dynamics_.l1NormImpulseDynamicsResidual(kkt_residual_);
  return violation;
}

} // namespace idocp