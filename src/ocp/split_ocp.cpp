#include "idocp/ocp/split_ocp.hpp"

#include <assert.h>


namespace idocp {

SplitOCP::SplitOCP(const Robot& robot, 
                   const std::shared_ptr<CostFunction>& cost,
                   const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(),
    kkt_residual_(robot),
    kkt_matrix_(robot),
    // robot_dynamics_(robot),
    contact_dynamics_(robot),
    riccati_gain_(robot),
    riccati_factorizer_(robot),
    s_tmp_(robot),
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
    riccati_gain_(),
    riccati_factorizer_(),
    s_tmp_(),
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


void SplitOCP::linearizeOCP(Robot& robot, const ContactStatus& contact_status,  
                            const double t, const double dtau, 
                            const Eigen::VectorXd& q_prev, 
                            const SplitSolution& s, 
                            const SplitSolution& s_next) {
  assert(dtau > 0);
  setContactStatusForKKT(contact_status);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  kkt_residual_.setZero();
  kkt_matrix_.setZero();
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  cost_->computeStageCostHessian(robot, cost_data_, t, dtau, s, kkt_matrix_);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s,
                                    kkt_residual_);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dtau, s, 
                                     kkt_matrix_, kkt_residual_);
  stateequation::LinearizeForwardEuler(robot, dtau, q_prev, s, s_next, 
                                       kkt_matrix_, kkt_residual_);
  contact_dynamics_.condenseContactDynamics(robot, contact_status, dtau, s, 
                                            kkt_matrix_, kkt_residual_);
}


void SplitOCP::backwardRiccatiRecursion(const double dtau,
                                        const RiccatiSolution& riccati_next, 
                                        RiccatiSolution& riccati) {
  assert(dtau > 0);
  riccati_factorizer_.factorizeBackwardRicursion(riccati_next, dtau, 
                                                 kkt_matrix_, kkt_residual_, 
                                                 riccati_gain_, riccati);
}


void SplitOCP::forwardRiccatiRecursion(const double dtau, SplitDirection& d,   
                                       SplitDirection& d_next) const {
  assert(dtau > 0);
  riccati_factorizer_.computeControlInputDirection(riccati_gain_, d);
  riccati_factorizer_.factorizeForwardRicursion(kkt_matrix_, kkt_residual_, d,
                                                dtau, d_next);
}


void SplitOCP::computeCondensedPrimalDirection(Robot& robot, const double dtau, 
                                               const RiccatiSolution& riccati,
                                               const SplitSolution& s, 
                                               SplitDirection& d) {
  assert(dtau > 0);
  riccati_factorizer_.computeCostateDirection(riccati, d);
  contact_dynamics_.computeCondensedPrimalDirection(robot, d);
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, dtau, s, d);
}


void SplitOCP::computeCondensedDualDirection(Robot& robot, const double dtau, 
                                             const SplitDirection& d_next, 
                                             SplitDirection& d) {
  assert(dtau > 0);
  contact_dynamics_.computeCondensedDualDirection(robot, dtau, kkt_matrix_,
                                                  kkt_residual_, 
                                                  d_next.dgmm(), d);
}


double SplitOCP::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


double SplitOCP::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


std::pair<double, double> SplitOCP::costAndConstraintViolation(
    Robot& robot, const double t, const double dtau, const SplitSolution& s) {
  assert(dtau > 0);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, s);
  cost += constraints_->costSlackBarrier(constraints_data_);
  return std::make_pair(cost, constraintViolation(dtau));
}


std::pair<double, double> SplitOCP::costAndConstraintViolation(
    Robot& robot, const ContactStatus& contact_status, const double step_size, 
    const double t, const double dtau, const SplitSolution& s, 
    const SplitDirection& d, const SplitSolution& s_next, 
    const SplitDirection& d_next) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  setContactStatusForKKT(contact_status);
  s_tmp_.setContactStatus(contact_status);
  s_tmp_.integratePrimal(robot, s, step_size, d);
  if (use_kinematics_) {
    robot.updateKinematics(s_tmp_.q, s_tmp_.v, s_tmp_.a);
  }
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, s_tmp_);
  cost += constraints_->costSlackBarrier(constraints_data_, step_size);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, dtau, 
                                             s_tmp_);
  stateequation::ComputeForwardEulerResidual(robot, step_size, dtau, s_tmp_,  
                                             s_next.q, s_next.v, d_next.dq(), 
                                             d_next.dv(), kkt_residual_);
  contact_dynamics_.computeContactDynamicsResidual(robot, contact_status, dtau, 
                                                   s_tmp_);
  return std::make_pair(cost, constraintViolation(dtau));
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


void SplitOCP::getStateFeedbackGain(Eigen::MatrixXd& Kq, 
                                    Eigen::MatrixXd& Kv) const {
  assert(Kq.cols() == riccati_gain_.Kq().cols());
  assert(Kq.rows() == riccati_gain_.Kq().rows());
  assert(Kv.cols() == riccati_gain_.Kv().cols());
  assert(Kv.rows() == riccati_gain_.Kv().rows());
  if (fd_like_elimination_) {
    Kq = riccati_gain_.Kq();
    Kv = riccati_gain_.Kv();
  }
  else {
    // robot_dynamics_.getStateFeedbackGain(riccati_gain_.Kq(), riccati_gain_.Kv(), 
    //                                      Kq, Kv);
  }
}


void SplitOCP::computeKKTResidual(Robot& robot, 
                                  const ContactStatus& contact_status, 
                                  const double t, const double dtau, 
                                  const Eigen::VectorXd& q_prev, 
                                  const SplitSolution& s,
                                  const SplitSolution& s_next) {
  assert(dtau > 0);
  setContactStatusForKKT(contact_status);
  kkt_residual_.setZero();
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, dtau, s);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s,
                                    kkt_residual_);
  stateequation::LinearizeForwardEuler(robot, dtau, q_prev, s, s_next, 
                                       kkt_matrix_, kkt_residual_);
  contact_dynamics_.linearizeContactDynamics(robot, contact_status, dtau, s, 
                                             kkt_matrix_, kkt_residual_);
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


double SplitOCP::constraintViolation(const double dtau) const {
  double violation = 0;
  violation += stateequation::L1NormStateEuqationResidual(kkt_residual_);
  violation += contact_dynamics_.l1NormContactDynamicsResidual(dtau);
  violation += constraints_->l1NormPrimalResidual(constraints_data_);
  return violation;
}

} // namespace idocp