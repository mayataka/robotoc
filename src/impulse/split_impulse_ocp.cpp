#include "idocp/impulse/split_impulse_ocp.hpp"

#include <assert.h>


namespace idocp {

ImpulseSplitOCP::ImpulseSplitOCP(
    const Robot& robot, const std::shared_ptr<ImpulseCostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(),
    kkt_residual_(robot),
    kkt_matrix_(robot),
    impulse_dynamics_(robot),
    riccati_factorizer_(robot),
    riccati_inverter_(robot),
    s_tmp_(robot),
    dimv_(robot.dimv()),
    dimf_(0),
    dimc_(0) {
}


ImpulseSplitOCP::ImpulseSplitOCP() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    kkt_residual_(),
    kkt_matrix_(),
    impulse_dynamics_(),
    riccati_factorizer_(),
    riccati_inverter_(),
    s_tmp_(),
    dimv_(0),
    dimf_(0),
    dimc_(0) {
}


ImpulseSplitOCP::~ImpulseSplitOCP() {
}


bool ImpulseSplitOCP::isFeasible(Robot& robot, const ImpulseSplitSolution& s) {
  return true;
  // return constraints_->isFeasible(robot, constraints_data_, s);
}


void ImpulseSplitOCP::initConstraints(Robot& robot,
                                      const ImpulseSplitSolution& s) { 
  // constraints_->setSlackAndDual(robot, constraints_data_, dtau, s);
}


void ImpulseSplitOCP::linearizeOCP(Robot& robot, 
                                   const ContactStatus& contact_status,  
                                   const double t,  
                                   const Eigen::VectorXd& q_prev, 
                                   const ImpulseSplitSolution& s, 
                                   const SplitSolution& s_next) {
  setContactStatusForKKT(contact_status);
  robot.updateKinematics(s.q, s.v);
  // condensing the impulse dynamics
  kkt_matrix.setZero();
  kkt_residual.setZero();
  cost_->computeStageCostHessian(robot, cost_data_, t, s, kkt_matrix_);
  cost_->computeStageCostDerivatives(robot, cost_data_, t, s, kkt_residual_);
  stateequation::LinearizeImpulseForwardEuler(robot, q_prev, s, s_next, 
                                              kkt_matrix_, kkt_residual_);
  impulse_dynamics_.condenseImpulseDynamics(robot, contact_status, s, 
                                            kkt_matrix_, kkt_residual_);
  kkt_matrix_.Qvq() = kkt_matrix_.Qqv().transpose();
}


void ImpulseSplitOCP::backwardRiccatiRecursion(
    const RiccatiFactorization& riccati_next, RiccatiFactorization& riccati) {
  riccati_factorizer_.factorize(kkt_matrix_, kkt_residual_, riccati_next, 
                                riccati);
}


void ImpulseSplitOCP::forwardRiccatiRecursion(ImpulseSplitDirection& d,   
                                              SplitDirection& d_next) {
  d_next.dq().noalias() = kkt_matrix.Fqq * d.dq() + kkt_residual_.Fq();
  d_next.dv().noalias() = kkt_matrix.Fvq * d.dq() 
                          + kkt_matrix.Fvv * d.dv() + kkt_residual_.Fv();
}


void ImpulseSplitOCP::computeCondensedDirection(Robot& robot, 
                                                const ImpulseSplitSolution& s, 
                                                const SplitDirection& d_next, 
                                                ImpulseSplitDirection& d) {
  impulse_dynamics_.computeCondensedDirection(kkt_matrix_, kkt_residual_, 
                                              d_next, d);
}

 
double ImpulseSplitOCP::maxPrimalStepSize() {
  return 1;
  // return constraints_->maxSlackStepSize(constraints_data_);
}


double ImpulseSplitOCP::maxDualStepSize() {
  return 1;
  // return constraints_->maxDualStepSize(constraints_data_);
}


std::pair<double, double> ImpulseSplitOCP::costAndConstraintViolation(
    Robot& robot, const double t, const ImpulseSplitSolution& s) {
  return std::make_pair(cost(robot, t, s), constraintViolation());
}


std::pair<double, double> ImpulseSplitOCP::costAndConstraintViolation(
    Robot& robot, const ContactStatus& contact_status, const double step_size, 
    const double t, const ImpulseSplitSolution& s, 
    const ImpulseSplitDirection& d, const SplitSolution& s_next, 
    const SplitDirection& d_next) {
  assert(step_size > 0);
  assert(step_size <= 1);
  setContactStatusForKKT(contact_status);
  s_tmp_.setContactStatus(contact_status);
  s_tmp_.dv = s.dv + step_size * d.ddv;
  s_tmp_.f_stack() = s.f_stack() + step_size * d.df();
  s_tmp_.set_f();
  robot.setContactForces(contact_status, s_tmp_.f);
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp_.q);
  s_tmp_.v = s.v + step_size * d.dv();
  robot.updateKinematics(s_tmp_.q, s_tmp_.v);
  // constraints_->computePrimalAndDualResidual(robot, constraints_data_, 
  //                                            s_tmp_);
  stateequation::ComputeImpulseForwardEulerResidual(robot, step_size, s_tmp_,  
                                                    s_next.q, s_next.v, 
                                                    d_next.dq(), d_next.dv(), 
                                                    kkt_residual_);
  impulse_dynamics_.computeImpulseDynamicsResidual(robot, contact_status, 
                                                   s_tmp_, kkt_residual_);
  return std::make_pair(cost(robot, t, s_tmp_), constraintViolation());
}


void ImpulseSplitOCP::updateDual(const double step_size) {
  assert(step_size > 0);
  assert(step_size <= 1);
  // constraints_->updateDual(constraints_data_, step_size);
}


void ImpulseSplitOCP::updatePrimal(Robot& robot, const double step_size, 
                                   const RiccatiFactorization& riccati, 
                                   const ImpulseSplitDirection& d, 
                                   ImpulseSplitSolution& s) {
  assert(step_size > 0);
  assert(step_size <= 1);
  s.lmd.noalias() 
      += step_size * (riccati.Pqq * d.dq() + riccati.Pqv * d.dv() - riccati.sq);
  s.gmm.noalias() 
      += step_size * (riccati.Pvq * d.dq() + riccati.Pvv * d.dv() - riccati.sv);
  robot.integrateConfiguration(d.dq(), step_size, s.q);
  s.v.noalias() += step_size * d.dv();
  s.dv.noalias() += step_size * d.ddv;
  s.f_stack().noalias() += step_size * d.df();
  s.set_f();
  s.beta.noalias() += step_size * d.dbeta;
  s.mu_stack().noalias() += step_size * d.dmu();
  s.set_mu_contact();
  // constraints_->updateSlack(constraints_data_, step_size);
}


void ImpulseSplitOCP::computeKKTResidual(Robot& robot, 
                                         const ContactStatus& contact_status, 
                                         const double t, 
                                         const Eigen::VectorXd& q_prev, 
                                         const SplitSolution& s,
                                         const SplitSolution& s_next) {
  setContactStatusForKKT(contact_status);
  kkt_residual_.setZero();
  robot.updateKinematics(s.q, s.v);
  cost_->computeStageCostDerivatives(robot, cost_data_, t, s, kkt_residual_);
  // constraints_->computePrimalAndDualResidual(robot, constraints_data_, dtau, s);
  // constraints_->augmentDualResidual(robot, constraints_data_, dtau, s,
  //                                   kkt_residual_);
  // constraints_->augmentDualResidual(robot, constraints_data_, dtau, s.u,
  //                                   kkt_residual_.lu);
  stateequation::LinearizeImpulseForwardEuler(robot, q_prev, s, s_next, 
                                              kkt_matrix_, kkt_residual_);
  impulse_dynamics_.linearizeImpulseDynamics(robot, contact_status, s, 
                                             kkt_matrix_, kkt_residual_);
}


double ImpulseSplitOCP::squaredNormKKTResidual(const double dtau) const {
  double error = 0;
  error += kkt_residual_.lq().squaredNorm();
  error += kkt_residual_.lv().squaredNorm();
  error += kkt_residual_.ldv.squaredNorm();
  error += kkt_residual_.lf().squaredNorm();
  error += stateequation::SquaredNormStateEuqationResidual(kkt_residual_);
  error += impulse_dynamics_.squaredNormImpulseDynamicsResidual(kkt_residual_);
  // error += constraints_->squaredNormPrimalAndDualResidual(constraints_data_);
  return error;
}


double ImpulseSplitOCP::cost(Robot& robot, const double t, 
                             const ImpulseSplitSolution& s) {
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, s);
  // cost += constraints_->costSlackBarrier(constraints_data_);
  return cost;
}


double ImpulseSplitOCP::constraintViolation() const {
  double violation = 0;
  violation += stateequation::L1NormStateEuqationResidual(kkt_residual_);
  violation += impulse_dynamics_.l1NormImpulseDynamicsResidual(kkt_residual_);
  // violation += constraints_->l1NormPrimalResidual(constraints_data_);
  return violation;
}

} // namespace idocp