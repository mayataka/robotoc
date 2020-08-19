#include "idocp/ocp/terminal_parnmpc.hpp"
#include "Eigen/LU"

#include <assert.h>


namespace idocp {

TerminalParNMPC::TerminalParNMPC(const Robot& robot, 
                                 const std::shared_ptr<CostFunction>& cost, 
                                 const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(robot),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot)),
    kkt_residual_(robot),
    kkt_matrix_(robot),
    kkt_matrix_inverse_(robot),
    state_equation_(robot),
    inverse_dynamics_(robot),
    x_res_(Eigen::VectorXd::Zero(2*robot.dimv())),
    dx_(Eigen::VectorXd::Zero(2*robot.dimv())) {
}


TerminalParNMPC::TerminalParNMPC() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    kkt_residual_(),
    kkt_matrix_(),
    kkt_matrix_inverse_(),
    state_equation_(),
    inverse_dynamics_(),
    x_res_(),
    dx_() {
}


TerminalParNMPC::~TerminalParNMPC() {
}


bool TerminalParNMPC::isFeasible(const Robot& robot, const SplitSolution& s) {
  return constraints_->isFeasible(robot, constraints_data_, s);
}


void TerminalParNMPC::initConstraints(const Robot& robot, const int time_step, 
                                   const double dtau, const SplitSolution& s) {
  assert(time_step >= 0);
  assert(dtau > 0);
  constraints_->setSlackAndDual(robot, constraints_data_, dtau, s);
}


void TerminalParNMPC::coarseUpdate(Robot& robot, const double t, 
                                   const double dtau, 
                                   const Eigen::VectorXd& q_prev, 
                                   const Eigen::VectorXd& v_prev,
                                   const SplitSolution& s, SplitDirection& d, 
                                   SplitSolution& s_new_coarse) {
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  assert(s.dimf() == robot.dimf());
  assert(s.dimc() == robot.dim_passive()+robot.dimf());
  assert(d.dimf() == robot.dimf());
  assert(d.dimc() == robot.dim_passive()+robot.dimf());
  assert(s_new_coarse.dimf() == robot.dimf());
  assert(s_new_coarse.dimc() == robot.dim_passive()+robot.dimf());
  kkt_residual_.setZero();
  kkt_residual_.setContactStatus(robot);
  kkt_matrix_.setZero();
  kkt_matrix_.setContactStatus(robot);
  kkt_matrix_inverse_.setZero();
  kkt_matrix_inverse_.setContactStatus(robot);
  if (robot.has_active_contacts()) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  computeKKTResidual(robot, t, dtau, q_prev, v_prev, s);
  cost_->computeStageCostHessian(robot, cost_data_, t, dtau, s, kkt_matrix_);
  cost_->computeTerminalCostHessian(robot, cost_data_, t, s, kkt_matrix_);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dtau, s, 
                                     kkt_matrix_, kkt_residual_);
  inverse_dynamics_.condenseInverseDynamics(kkt_matrix_, kkt_residual_);
  inverse_dynamics_.condenseEqualityConstraint(dtau, kkt_matrix_, kkt_residual_);
  // coarse update of the solution
  kkt_matrix_.symmetrize();
  kkt_matrix_.invert(kkt_matrix_inverse_.KKT_matrix_inverse());
  d.split_direction() 
      = kkt_matrix_inverse_.KKT_matrix_inverse() * kkt_residual_.KKT_residual();
  s_new_coarse.lmd = s.lmd - d.dlmd();
  s_new_coarse.gmm = s.gmm - d.dgmm();
  s_new_coarse.mu_active() = s.mu_active() - d.dmu();
  s_new_coarse.a = s.a - d.da();
  s_new_coarse.f_active() = s.f_active() - d.df();
  robot.integrateConfiguration(s.q, d.dq(), -1, s_new_coarse.q);
  s_new_coarse.v = s.v - d.dv();
}


void TerminalParNMPC::getAuxiliaryMatrix(Eigen::MatrixXd& aux_mat) const {
  aux_mat = - kkt_matrix_inverse_.auxiliaryMatrix();
}


void TerminalParNMPC::backwardCorrectionSerial(const Robot& robot,
                                               const SplitSolution& s_old_next, 
                                               const SplitSolution& s_new_next,
                                               SplitSolution& s_new) {
  x_res_.head(robot.dimv()) = s_new_next.lmd - s_old_next.lmd;
  x_res_.tail(robot.dimv()) = s_new_next.gmm - s_old_next.gmm;
  dx_ = kkt_matrix_inverse_.backwardCorrectionSerialCoeff() * x_res_;
  s_new.lmd.noalias() -= dx_.head(robot.dimv());
  s_new.gmm.noalias() -= dx_.tail(robot.dimv());
}


void TerminalParNMPC::backwardCorrectionParallel(const Robot& robot, 
                                                 SplitDirection& d,
                                                 SplitSolution& s_new) {
  d.backwardCorrectionParallelDirection() 
      = kkt_matrix_inverse_.backwardCorrectionParallelCoeff() * x_res_;
  s_new.mu_active().noalias() -= d.dmu();
  s_new.a.noalias() -= d.da();
  s_new.f_active().noalias() -= d.df();
  robot.integrateConfiguration(d.dq(), -1, s_new.q);
  s_new.v.noalias() -= d.dv();
}


void TerminalParNMPC::forwardCorrectionSerial(const Robot& robot,
                                              const SplitSolution& s_old_prev,
                                              const SplitSolution& s_new_prev, 
                                              SplitSolution& s_new) {
  robot.subtractConfiguration(s_new_prev.q, s_old_prev.q, 
                              x_res_.head(robot.dimv()));
  x_res_.tail(robot.dimv()) = s_new_prev.v - s_old_prev.v;
  dx_ = kkt_matrix_inverse_.forwardCorrectionSerialCoeff() * x_res_;
  robot.integrateConfiguration(dx_.head(robot.dimv()), -1, s_new.q);
  s_new.v.noalias() -= dx_.tail(robot.dimv());
}


void TerminalParNMPC::forwardCorrectionParallel(const Robot& robot,
                                                SplitDirection& d,
                                                SplitSolution& s_new) {
  d.forwardCorrectionParallelDirection() 
      = kkt_matrix_inverse_.forwardCorrectionParallelCoeff() * x_res_;
  s_new.lmd.noalias() -= d.dlmd();
  s_new.gmm.noalias() -= d.dgmm();
  s_new.mu_active().noalias() -= d.dmu();
  s_new.a.noalias() -= d.da();
  s_new.f_active().noalias() -= d.df();
}


void TerminalParNMPC::computePrimalAndDualDirection(const Robot& robot, 
                                                    const double dtau, 
                                                    const SplitSolution& s, 
                                                    const SplitSolution& s_new,
                                                    SplitDirection& d) {
  d.dlmd() = s_new.lmd - s.lmd;
  d.dgmm() = s_new.gmm - s.gmm;
  d.dmu() = s_new.mu_active() - s.mu_active();
  d.da() = s_new.a - s.a;
  d.df() = s_new.f_active() - s.f_active();
  robot.subtractConfiguration(s_new.q, s.q, d.dq());
  d.dv() = s_new.v - s.v;
  inverse_dynamics_.computeCondensedDirection(dtau, kkt_matrix_, kkt_residual_, d);
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, dtau, d);
}

 
double TerminalParNMPC::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


double TerminalParNMPC::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


std::pair<double, double> TerminalParNMPC::costAndConstraintsViolation(
    Robot& robot, const double t, const double dtau, const SplitSolution& s) {
  assert(dtau > 0);
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, s);
  cost += cost_->phi(robot, cost_data_, t, s);
  cost += constraints_->costSlackBarrier(constraints_data_);
  double violation = 0;
  violation += state_equation_.violationL1Norm(kkt_residual_);
  violation += inverse_dynamics_.violationL1Norm(dtau, kkt_residual_);
  violation += equalityconstraints::ViolationL1Norm(kkt_residual_);
  violation += constraints_->residualL1Nrom(robot, constraints_data_, dtau, s);
  return std::make_pair(cost, violation);
}


std::pair<double, double> TerminalParNMPC::costAndConstraintsViolation(
    Robot& robot, const double step_size, const double t, const double dtau, 
    const SplitSolution& s_prev, const SplitDirection& d_prev,
    const SplitSolution& s, const SplitDirection& d, SplitSolution& s_tmp) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp.q);
  s_tmp.v = s.v + step_size * d.dv();
  s_tmp.a = s.a + step_size * d.da();
  if (robot.has_active_contacts()) {
    s_tmp.f_active() = s.f_active() + step_size * d.df();
    robot.setContactForces(s_tmp.f);
  }
  s_tmp.u = s.u + step_size * d.du;
  if (robot.has_active_contacts()) {
    robot.updateKinematics(s_tmp.q, s_tmp.v, s_tmp.a);
  }
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, s_tmp);
  cost += cost_->phi(robot, cost_data_, t, s_tmp);
  cost += constraints_->costSlackBarrier(constraints_data_, step_size);
  double violation = 0;
  violation += state_equation_.violationL1Norm(robot, step_size, dtau, s_prev.q, 
                                               s_prev.v, d_prev.dq(), d_prev.dv(),
                                               s, kkt_residual_);
  violation += inverse_dynamics_.violationL1Norm(robot, dtau, s, kkt_residual_);
  violation += equalityconstraints::ViolationL1Norm(robot, dtau, s, kkt_residual_);
  violation += constraints_->residualL1Nrom(robot, constraints_data_, dtau, s);
  return std::make_pair(cost, violation);
}


void TerminalParNMPC::updateDual(const double step_size) {
  assert(step_size > 0);
  assert(step_size <= 1);
  constraints_->updateDual(constraints_data_, step_size);
}


void TerminalParNMPC::updatePrimal(Robot& robot, const double step_size, 
                                const double dtau, const SplitDirection& d, 
                                SplitSolution& s) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  s.lmd.noalias() += step_size * d.dlmd();
  s.gmm.noalias() += step_size * d.dgmm();
  s.mu_active().noalias() += step_size * d.dmu();
  s.a.noalias() += step_size * d.da();
  s.f_active().noalias() += step_size * d.df();
  robot.integrateConfiguration(d.dq(), step_size, s.q);
  s.v.noalias() += step_size * d.dv();
  s.u.noalias() += step_size * d.du;
  s.beta.noalias() += step_size * d.dbeta;
  constraints_->updateSlack(constraints_data_, step_size);
}


void TerminalParNMPC::getStateFeedbackGain(Eigen::MatrixXd& Kq, 
                                        Eigen::MatrixXd& Kv) const {
  // Kq = du_dq_ + du_da_ * Kaq_ + du_df_.leftCols(dimf_) * Kfq_.topRows(dimf_);
  // Kv = du_dv_ + du_da_ * Kav_ + du_df_.leftCols(dimf_) * Kfv_.topRows(dimf_);
}


double TerminalParNMPC::squaredKKTErrorNorm(Robot& robot, const double t, 
                                            const double dtau, 
                                            const Eigen::VectorXd& q_prev, 
                                            const Eigen::VectorXd& v_prev, 
                                            const SplitSolution& s) {
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  kkt_matrix_.setZero();
  kkt_matrix_.setContactStatus(robot);
  kkt_residual_.setZero();
  kkt_residual_.setContactStatus(robot);
  if (robot.has_active_contacts()) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  computeKKTResidual(robot, t, dtau, q_prev, v_prev, s);
  double error = kkt_residual_.squaredKKTErrorNorm(dtau);
  error += constraints_->squaredKKTErrorNorm(robot, constraints_data_, dtau, s);
  return error;
}


void TerminalParNMPC::computeKKTResidual(Robot& robot, const double t, 
                                         const double dtau, 
                                         const Eigen::VectorXd& q_prev, 
                                         const Eigen::VectorXd& v_prev, 
                                         const SplitSolution& s) {
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  cost_->computeTerminalCostDerivatives(robot, cost_data_, t, s, 
                                        kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, 
                                    kkt_residual_);
  state_equation_.linearizeStateEquationTerminal(robot, dtau, q_prev, v_prev, s, 
                                                 kkt_matrix_, kkt_residual_);
  equalityconstraints::LinearizeEqualityConstraints(robot, dtau, s, 
                                                    kkt_matrix_, kkt_residual_);
  inverse_dynamics_.linearizeInverseDynamics(robot, dtau, s, kkt_residual_);
}


} // namespace idocp