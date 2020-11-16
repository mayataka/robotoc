#include "idocp/ocp/split_ocp.hpp"

#include <cassert>


namespace idocp {

SplitOCP::SplitOCP(const Robot& robot, 
                   const std::shared_ptr<CostFunction>& cost,
                   const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot)),
    kkt_residual_(robot),
    kkt_matrix_(robot),
    // robot_dynamics_(robot),
    contact_dynamics_(robot),
    riccati_factorizer_(robot),
    has_floating_base_(robot.has_floating_base()),
    use_kinematics_(false),
    fd_like_elimination_(false) {
  if (cost_->useKinematics() || constraints_->useKinematics() 
                             || robot.max_point_contacts() > 0) {
    use_kinematics_ = true;
  }
  if (robot.has_floating_base()) {
    fd_like_elimination_ = true;
  }
}


SplitOCP::SplitOCP() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    kkt_residual_(),
    kkt_matrix_(),
    // robot_dynamics_(),
    contact_dynamics_(),
    riccati_factorizer_(),
    has_floating_base_(false),
    use_kinematics_(false),
    fd_like_elimination_(false) {
}


SplitOCP::~SplitOCP() {
}


bool SplitOCP::isFeasible(Robot& robot, const SplitSolution& s) {
  return constraints_->isFeasible(robot, constraints_data_, s);
}


void SplitOCP::initConstraints(Robot& robot, const int time_step, 
                               const double dtau, const SplitSolution& s) { 
  assert(time_step >= 0);
  assert(dtau > 0);
  constraints_data_ = constraints_->createConstraintsData(robot, time_step);
  constraints_->setSlackAndDual(robot, constraints_data_, dtau, s);
}


void SplitOCP::backwardRiccatiRecursion(
    const RiccatiFactorization& riccati_next, const double dtau, 
    RiccatiFactorization& riccati) {
  assert(dtau > 0);
  riccati_factorizer_.backwardRiccatiRecursion(riccati_next, dtau, kkt_matrix_, 
                                               kkt_residual_, riccati);
}


void SplitOCP::forwardRiccatiRecursionParallel(
    const bool exist_state_constraint) {
  riccati_factorizer_.forwardRiccatiRecursionParallel(kkt_matrix_, 
                                                      kkt_residual_,
                                                      exist_state_constraint);
}


void SplitOCP::forwardRiccatiRecursionSerial(
    const RiccatiFactorization& riccati, const double dtau,
    RiccatiFactorization& riccati_next, const bool exist_state_constraint) {
  assert(dtau > 0);
  riccati_factorizer_.forwardRiccatiRecursionSerial(riccati, kkt_matrix_, 
                                                    kkt_residual_, dtau,
                                                    riccati_next,
                                                    exist_state_constraint);
}


void SplitOCP::forwardRiccatiRecursionUnconstrained(
    SplitDirection& d, const double dtau, SplitDirection& d_next) const {
  riccati_factorizer_.forwardRiccatiRecursionUnconstrained(kkt_matrix_, 
                                                           kkt_residual_, d, 
                                                           dtau, d_next);
}


void SplitOCP::computePrimalDirection(
    Robot& robot, const double dtau, const RiccatiFactorization& riccati, 
    const SplitSolution& s, SplitDirection& d, 
    const bool exist_state_constraint) {
  riccati_factorizer_.computeCostateDirection(riccati, d, 
                                              exist_state_constraint);
  riccati_factorizer_.computeControlInputDirection(riccati, d, 
                                                   exist_state_constraint);
  contact_dynamics_.computeCondensedPrimalDirection(robot, d);
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, dtau, 
                                             s, d);
}

double SplitOCP::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


double SplitOCP::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


void SplitOCP::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(constraints_data_, dual_step_size);
}


void SplitOCP::updatePrimal(Robot& robot, const double primal_step_size, 
                            const double dtau, const SplitDirection& d, 
                            SplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  assert(dtau > 0);
  s.integrate(robot, primal_step_size, d);
  constraints_->updateSlack(constraints_data_, primal_step_size);
}


double SplitOCP::squaredNormKKTResidual(const double dtau) const {
  double error = 0;
  error += kkt_residual_.lx().squaredNorm();
  error += kkt_residual_.la.squaredNorm();
  error += kkt_residual_.lf().squaredNorm();
  if (has_floating_base_) {
    error += kkt_residual_.lu_passive.squaredNorm();
  }
  error += kkt_residual_.lu().squaredNorm();
  error += stateequation::SquaredNormStateEuqationResidual(kkt_residual_);
  error += contact_dynamics_.squaredNormContactDynamicsResidual(dtau);
  error += constraints_->squaredNormPrimalAndDualResidual(constraints_data_);
  return error;
}


double SplitOCP::stageCost(Robot& robot, const double t, const double dtau, 
                           const SplitSolution& s, 
                           const double primal_step_size) {
  assert(dtau > 0);
  assert(primal_step_size >= 0);
  assert(primal_step_size <= 1);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, s);
  if (primal_step_size > 0) {
    cost += constraints_->costSlackBarrier(constraints_data_, primal_step_size);
  }
  else {
    cost += constraints_->costSlackBarrier(constraints_data_);
  }
  return cost;
}


double SplitOCP::constraintViolation(Robot& robot, 
                                     const ContactStatus& contact_status, 
                                     const double t, const double dtau, 
                                     const SplitSolution& s, 
                                     const Eigen::VectorXd& q_next, 
                                     const Eigen::VectorXd& v_next) {
  setContactStatusForKKT(contact_status);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, dtau, s);
  stateequation::ComputeForwardEulerResidual(robot, dtau, s, q_next, v_next, 
                                             kkt_residual_);
  contact_dynamics_.computeContactDynamicsResidual(robot, contact_status, dtau, s);
  double violation = 0;
  violation += constraints_->l1NormPrimalResidual(constraints_data_);
  violation += stateequation::L1NormStateEuqationResidual(kkt_residual_);
  violation += contact_dynamics_.l1NormContactDynamicsResidual(dtau);
  return violation;
}

} // namespace idocp