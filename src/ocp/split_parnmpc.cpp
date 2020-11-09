#include "idocp/ocp/split_parnmpc.hpp"

#include <cassert>


namespace idocp {

SplitParNMPC::SplitParNMPC(const Robot& robot, 
                           const std::shared_ptr<CostFunction>& cost, 
                           const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(robot),
    constraints_(constraints),
    constraints_data_(),
    kkt_residual_(robot),
    kkt_matrix_(robot),
    contact_dynamics_(robot),
    backward_correction_(robot),
    use_kinematics_(false),
    has_floating_base_(robot.has_floating_base()) {
  if (cost_->useKinematics() || constraints_->useKinematics() 
                             || robot.max_point_contacts() > 0) {
    use_kinematics_ = true;
  }
}


SplitParNMPC::SplitParNMPC() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    kkt_residual_(),
    kkt_matrix_(),
    contact_dynamics_(),
    backward_correction_(),
    use_kinematics_(false) {
}


SplitParNMPC::~SplitParNMPC() {
}


bool SplitParNMPC::isFeasible(Robot& robot, const SplitSolution& s) {
  return constraints_->isFeasible(robot, constraints_data_, s);
}


void SplitParNMPC::initConstraints(Robot& robot, const int time_step, 
                                   const double dtau, const SplitSolution& s) {
  assert(time_step >= 0);
  assert(dtau > 0);
  constraints_data_ = constraints_->createConstraintsData(robot, time_step);
  constraints_->setSlackAndDual(robot, constraints_data_, dtau, s);
}


void SplitParNMPC::linearizeOCP(Robot& robot, 
                                const ContactStatus& contact_status,
                                const double t, const double dtau, 
                                const Eigen::VectorXd& q_prev, 
                                const Eigen::VectorXd& v_prev,
                                const SplitSolution& s,
                                const SplitSolution& s_next) {
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  setContactStatusForKKT(contact_status);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  kkt_residual_.setZero();
  kkt_matrix_.setZero();
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s,
                                    kkt_residual_);
  stateequation::LinearizeBackwardEuler(robot, dtau, q_prev, v_prev, s, s_next, 
                                        kkt_matrix_, kkt_residual_);
  contact_dynamics_.linearizeContactDynamics(robot, contact_status, dtau, s, 
                                             kkt_matrix_, kkt_residual_);
  cost_->computeStageCostHessian(robot, cost_data_, t, dtau, s, kkt_matrix_);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dtau, s, 
                                     kkt_matrix_, kkt_residual_);
  contact_dynamics_.condenseContactDynamics(robot, contact_status, dtau, 
                                            kkt_matrix_, kkt_residual_);
}


void SplitParNMPC::coarseUpdate(const Robot& robot, const SplitSolution& s,
                                const Eigen::MatrixXd& aux_mat_next_old,
                                SplitDirection& d, 
                                SplitSolution& s_new_coarse) {
  kkt_matrix_.Qxx().noalias() += aux_mat_next_old;
  kkt_matrix_.symmetrize(); 
  backward_correction_.coarseUpdate(robot, s, d, kkt_matrix_, kkt_residual_, 
                                    s_new_coarse);
}


void SplitParNMPC::getAuxiliaryMatrix(Eigen::MatrixXd& aux_mat) const {
  backward_correction_.getAuxiliaryMatrix(aux_mat);
}


void SplitParNMPC::backwardCorrectionSerial(const SplitSolution& s_next, 
                                            const SplitSolution& s_new_next,
                                            SplitSolution& s_new) {
  backward_correction_.backwardCorrectionSerial(s_next, s_new_next, s_new);
}


void SplitParNMPC::backwardCorrectionParallel(const Robot& robot, 
                                              SplitDirection& d,
                                              SplitSolution& s_new) const {
  backward_correction_.backwardCorrectionParallel(robot, d, s_new);
}


void SplitParNMPC::forwardCorrectionSerial(const Robot& robot,
                                           const SplitSolution& s_prev,
                                           const SplitSolution& s_new_prev, 
                                           SplitSolution& s_new) {
  backward_correction_.forwardCorrectionSerial(robot, s_prev, s_new_prev, s_new);
}


void SplitParNMPC::forwardCorrectionParallel(SplitDirection& d,
                                             SplitSolution& s_new) const {
  backward_correction_.forwardCorrectionParallel(d, s_new);
}


void SplitParNMPC::computePrimalAndDualDirection(Robot& robot, 
                                                 const double dtau, 
                                                 const SplitSolution& s, 
                                                 const SplitSolution& s_new,
                                                 SplitDirection& d) {
  backward_correction_.computeDirection(robot, s, s_new, d);
  contact_dynamics_.computeCondensedPrimalDirection(robot, d);
  contact_dynamics_.computeCondensedDualDirection(robot, dtau, kkt_matrix_,
                                                  kkt_residual_, d.dgmm(), d);
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, dtau, s, d);
}

 
double SplitParNMPC::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


double SplitParNMPC::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


void SplitParNMPC::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(constraints_data_, dual_step_size);
}


void SplitParNMPC::updatePrimal(Robot& robot, const double primal_step_size, 
                                const double dtau, const SplitDirection& d, 
                                SplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  assert(dtau > 0);
  s.integrate(robot, primal_step_size, d);
  constraints_->updateSlack(constraints_data_, primal_step_size);
}


void SplitParNMPC::computeKKTResidual(Robot& robot, 
                                      const ContactStatus& contact_status, 
                                      const double t, const double dtau, 
                                      const Eigen::VectorXd& q_prev, 
                                      const Eigen::VectorXd& v_prev, 
                                      const SplitSolution& s, 
                                      const SplitSolution& s_next) {
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
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
  stateequation::LinearizeBackwardEuler(robot, dtau, q_prev, v_prev, s, s_next,
                                        kkt_matrix_, kkt_residual_);
  contact_dynamics_.linearizeContactDynamics(robot, contact_status, dtau, s, 
                                             kkt_matrix_, kkt_residual_);
}


double SplitParNMPC::squaredNormKKTResidual(const double dtau) const {
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


double SplitParNMPC::stageCost(Robot& robot, const double t, const double dtau, 
                               const SplitSolution& s,
                               const double primal_step_size) {
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


double SplitParNMPC::constraintViolation(Robot& robot, 
                                         const ContactStatus& contact_status, 
                                         const double t, const double dtau, 
                                         const Eigen::VectorXd& q_prev, 
                                         const Eigen::VectorXd& v_prev,
                                         const SplitSolution& s) {
  setContactStatusForKKT(contact_status);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, dtau, s);
  stateequation::ComputeBackwardEulerResidual(robot, dtau, q_prev, v_prev, s, 
                                              kkt_residual_);
  contact_dynamics_.computeContactDynamicsResidual(robot, contact_status, dtau, s);
  double violation = 0;
  violation += constraints_->l1NormPrimalResidual(constraints_data_);
  violation += stateequation::L1NormStateEuqationResidual(kkt_residual_);
  violation += contact_dynamics_.l1NormContactDynamicsResidual(dtau);
  return violation;
}


void SplitParNMPC::getStateFeedbackGain(Eigen::MatrixXd& Kq, 
                                        Eigen::MatrixXd& Kv) const {
  // if (fd_like_elimination_) {
  //   Kq = riccati_gain_.Kq();
  //   Kv = riccati_gain_.Kv();
  // }
  // else {
  //   // robot_dynamics_.getStateFeedbackGain(riccati_gain_.Kq(), riccati_gain_.Kv(), 
  //   //                                      Kq, Kv);
  // }
}

} // namespace idocp