#include "idocp/ocp/split_parnmpc.hpp"
#include "idocp/ocp/parnmpc_linearizer.hpp"
#include "Eigen/LU"

#include <assert.h>


namespace idocp {

SplitParNMPC::SplitParNMPC(const Robot& robot, 
                           const std::shared_ptr<CostFunction>& cost, 
                           const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(robot),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot)),
    kkt_residual_(robot),
    kkt_matrix_(robot),
    kkt_matrix_inverse_(robot),
    linearizer_(robot),
    id_condenser_(robot),
    x_res_(Eigen::VectorXd::Zero(2*robot.dimv())),
    dx_(Eigen::VectorXd::Zero(2*robot.dimv())) {
}


SplitParNMPC::SplitParNMPC() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    kkt_residual_(),
    kkt_matrix_(),
    kkt_matrix_inverse_(),
    linearizer_(),
    id_condenser_(),
    x_res_(),
    dx_() {
}


SplitParNMPC::~SplitParNMPC() {
}


bool SplitParNMPC::isFeasible(const Robot& robot, const SplitSolution& s) {
  return constraints_->isFeasible(robot, constraints_data_, 
                                  s.a, s.f, s.q, s.v, s.u);
}


void SplitParNMPC::initConstraints(const Robot& robot, const int time_step, 
                                   const double dtau, const SplitSolution& s) {
  assert(time_step >= 0);
  assert(dtau > 0);
  constraints_->setSlackAndDual(robot, constraints_data_, dtau, 
                                s.a, s.f, s.q, s.v, s.u);
}


void SplitParNMPC::coarseUpdate(Robot& robot, const double t, const double dtau, 
                                const Eigen::VectorXd& q_prev, 
                                const Eigen::VectorXd& v_prev,
                                const SplitSolution& s,
                                const Eigen::VectorXd& lmd_next,
                                const Eigen::VectorXd& gmm_next,
                                const Eigen::VectorXd& q_next,
                                const Eigen::MatrixXd& aux_mat_next_old,
                                SplitDirection& d, 
                                SplitSolution& s_new_coarse) {
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  assert(lmd_next.size() == robot.dimv());
  assert(gmm_next.size() == robot.dimv());
  assert(q_next.size() == robot.dimq());
  assert(aux_mat_next_old.rows() == 2*robot.dimv());
  assert(aux_mat_next_old.cols() == 2*robot.dimv());
  assert(s.dimc() == robot.dim_passive()+robot.dimf());
  assert(d.dimc() == robot.dim_passive()+robot.dimf());
  assert(s_new_coarse.dimc() == robot.dim_passive()+robot.dimf());
  assert(s.dimf() == robot.dimf());
  assert(d.dimf() == robot.dimf());
  assert(s_new_coarse.dimf() == robot.dimf());
  kkt_residual_.setZero();
  kkt_residual_.setContactStatus(robot);
  kkt_matrix_.setZero();
  kkt_matrix_.setContactStatus(robot);
  kkt_residual_inverse_.setZero();
  kkt_residual_inverse_.setContactStatus(robot);
  if (robot.dimf() > 0) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  // Forms the Newton residual by linearizing cost, dynamics, and constraints. 
  linearizer_.linearizeCostAndConstraints(robot, cost_, cost_data_, 
                                          constraints_, constraints_data_, 
                                          t, dtau, s, kkt_residual_);
  linearizer_.linearizeStateEquation(robot, dtau, q_prev, v_prev, s, lmd_next, 
                                     gmm_next, q_next, kkt_residual_, kkt_matrix_);
  linearizer_.linearizeContactConstraints(robot, dtau, kkt_residual_, kkt_matrix_);
  // Newton residual with respect to inverse dynamics.
  id_condenser_.setContactStatus(robot);
  id_condenser_.linearizeCostAndConstraints(robot, cost_, cost_data_, 
                                            constraints_, constraints_data_, 
                                            t, dtau, s);
  id_condenser_.condenseInequalityConstraints(robot, constraints_, 
                                              constraints_data_, t, dtau, s);
  id_condenser_.linearizeInverseDynamics(robot, s);
  id_condenser_.augmentInverseDynamicsDerivatives(dtau, s, kkt_residual_);
  id_condenser_.condenseInverseDynamics(kkt_residual_, kkt_matrix_);
  if (robot.has_floating_base()) {
    id_condenser_.condenseFloatingBaseConstraint(dtau, s, kkt_residual_, kkt_matrix_);
  }
  // Augment the equality constraints 
  kkt_residual_.lq().noalias() += kkt_matrix_.Cq().transpose() * s.mu_active();
  kkt_residual_.lv().noalias() += kkt_matrix_.Cv().transpose() * s.mu_active();
  kkt_residual_.la().noalias() += kkt_matrix_.Ca().transpose() * s.mu_active();
  if (robot.dimf() > 0) {
    kkt_residual_.lf().noalias() += kkt_matrix_.Cf().transpose() * s.mu_active();
  }
  // Condense the slack and dual variables of the inequality constraints 
  constraints_->condenseSlackAndDual(robot, constraints_data_, dtau, 
                                     s.a, s.f, s.q, s.v, kkt_matrix_.Qaa(),
                                     kkt_matrix_.Qff(), kkt_matrix_.Qqq(),
                                     kkt_matrix_.Qvv(), kkt_residual_.la(), 
                                     kkt_residual_.lf(), kkt_residual_.lq(), 
                                     kkt_residual_.lv());
  // Augment the cost function Hessian. 
  cost_->augment_lqq(robot, cost_data_, t, dtau, s.q, s.v, s.a, kkt_matrix_.Qqq());
  cost_->augment_lvv(robot, cost_data_, t, dtau, s.q, s.v, s.a, kkt_matrix_.Qvv());
  cost_->augment_laa(robot, cost_data_, t, dtau, s.q, s.v, s.a, kkt_matrix_.Qaa());
  if (robot.dimf() > 0) {
    cost_->augment_lff(robot, cost_data_, t, dtau, s.f, kkt_matrix_.Qff());
  }
  kkt_matrix_.Qxx().noalias() += aux_mat_next_old;
  kkt_matrix_.symmetrize();
  const int dim_kkt = kkt_matrix_.dimKKT();
  kkt_matrix_.invert(kkt_matrix_inverse_.KKT_matrix_inverse());
  // coarse update of the solution
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


void SplitParNMPC::coarseUpdate(Robot& robot, const double t, const double dtau, 
                                const Eigen::VectorXd& q_prev, 
                                const Eigen::VectorXd& v_prev,
                                const SplitSolution& s, SplitDirection& d, 
                                SplitSolution& s_new_coarse) {
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  assert(s.dimc() == robot.dim_passive()+robot.dimf());
  assert(d.dimc() == robot.dim_passive()+robot.dimf());
  assert(s_new_coarse.dimc() == robot.dim_passive()+robot.dimf());
  assert(s.dimf() == robot.dimf());
  assert(d.dimf() == robot.dimf());
  assert(s_new_coarse.dimf() == robot.dimf());
  kkt_residual_.setZero();
  kkt_residual_.setContactStatus(robot);
  kkt_matrix_.setZero();
  kkt_matrix_.setContactStatus(robot);
  kkt_residual_inverse_.setZero();
  kkt_residual_inverse_.setContactStatus(robot);
  if (robot.dimf() > 0) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  // Forms the Newton residual by linearizing cost, dynamics, and constraints. 
  linearizer_.linearizeCostAndConstraints(robot, cost_, cost_data_, 
                                          constraints_, constraints_data_, 
                                          t, dtau, s, kkt_residual_);
  linearizer_.linearizeStateEquation(robot, dtau, q_prev, v_prev, s, 
                                     kkt_residual_, kkt_matrix_);
  linearizer_.linearizeContactConstraints(robot, dtau, kkt_residual_, kkt_matrix_);
  // Augmnet the partial derivatives of the inequality constriants.
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, 
                                    kkt_residual_.la(), 
                                    kkt_residual_.lf(), 
                                    kkt_residual_.lq(), 
                                    kkt_residual_.lv());
  // Newton residual with respect to inverse dynamics.
  id_condenser_.setContactStatus(robot);
  id_condenser_.linearizeCostAndConstraints(robot, cost_, cost_data_, 
                                            constraints_, constraints_data_, 
                                            t, dtau, s);
  id_condenser_.condenseInequalityConstraints(robot, constraints_, 
                                              constraints_data_, t, dtau, s);
  id_condenser_.linearizeInverseDynamics(robot, s);
  id_condenser_.augmentInverseDynamicsDerivatives(dtau, s, kkt_residual_);
  id_condenser_.condenseInverseDynamics(kkt_residual_, kkt_matrix_);
  if (robot.has_floating_base()) {
    id_condenser_.condenseFloatingBaseConstraint(dtau, s, kkt_residual_, kkt_matrix_);
  }
  // Augment the equality constraints 
  kkt_residual_.lq().noalias() += kkt_matrix_.Cq().transpose() * s.mu_active();
  kkt_residual_.lv().noalias() += kkt_matrix_.Cv().transpose() * s.mu_active();
  kkt_residual_.la().noalias() += kkt_matrix_.Ca().transpose() * s.mu_active();
  if (robot.dimf() > 0) {
    kkt_residual_.lf().noalias() += kkt_matrix_.Cf().transpose() * s.mu_active();
  }
  // Condense the slack and dual variables of the inequality constraints 
  constraints_->condenseSlackAndDual(robot, constraints_data_, dtau, 
                                     s.a, s.f, s.q, s.v, kkt_matrix_.Qaa(),
                                     kkt_matrix_.Qff(), kkt_matrix_.Qqq(),
                                     kkt_matrix_.Qvv(), kkt_residual_.la(), 
                                     kkt_residual_.lf(), kkt_residual_.lq(), 
                                     kkt_residual_.lv());
  // Augment the cost function Hessian. 
  cost_->augment_lqq(robot, cost_data_, t, dtau, s.q, s.v, s.a, kkt_matrix_.Qqq());
  cost_->augment_lvv(robot, cost_data_, t, dtau, s.q, s.v, s.a, kkt_matrix_.Qvv());
  cost_->augment_laa(robot, cost_data_, t, dtau, s.q, s.v, s.a, kkt_matrix_.Qaa());
  if (robot.dimf() > 0) {
    cost_->augment_lff(robot, cost_data_, t, dtau, s.f, kkt_matrix_.Qff());
  }
  cost_->augment_phiqq(robot, cost_data_, t, s.q, s.v, kkt_matrix_.Qqq());
  cost_->augment_phivv(robot, cost_data_, t, s.q, s.v, kkt_matrix_.Qvv());
  kkt_matrix_.symmetrize();
  std::cout << kkt_matrix_.KKT_matrix() << std::endl;
  const int dim_kkt = kkt_matrix_.dimKKT();
  kkt_matrix_.invert(kkt_matrix_inverse_.KKT_matrix_inverse());
  // coarse update of the solution
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


void SplitParNMPC::getAuxiliaryMatrix(Eigen::MatrixXd& aux_mat) const {
  aux_mat = - kkt_matrix_inverse_.auxiliaryMatrix();
}


void SplitParNMPC::backwardCorrectionSerial(const Robot& robot,
                                            const SplitSolution& s_old_next, 
                                            const SplitSolution& s_new_next,
                                            SplitSolution& s_new) {
  x_res_.head(robot.dimv()) = s_new_next.lmd - s_old_next.lmd;
  x_res_.tail(robot.dimv()) = s_new_next.gmm - s_old_next.gmm;
  dx_ = kkt_matrix_inverse_.backwardCorrectionSerialCoeff() * x_res_;
  s_new.lmd.noalias() -= dx_.head(robot.dimv());
  s_new.gmm.noalias() -= dx_.tail(robot.dimv());
}


void SplitParNMPC::backwardCorrectionParallel(const Robot& robot, 
                                              SplitDirection& d,
                                              SplitSolution& s_new) {
  const int dim_kkt = d.dimKKT();
  const int dimx = 2*robot.dimv();
  d.split_direction().segment(dimx, dim_kkt-dimx) = kkt_matrix_inverse_.backwardCorrectionParallelCoeff() * x_res_;
  s_new.mu_active().noalias() -= d.dmu();
  s_new.a.noalias() -= d.da();
  s_new.f_active().noalias() -= d.df();
  robot.integrateConfiguration(d.dq(), -1, s_new.q);
  s_new.v.noalias() -= d.dv();
}


void SplitParNMPC::forwardCorrectionSerial(const Robot& robot,
                                           const SplitSolution& s_old_prev,
                                           const SplitSolution& s_new_prev, 
                                           SplitSolution& s_new) {
  robot.subtractConfiguration(s_new_prev.q, s_old_prev.q, x_res_.head(robot.dimv()));
  x_res_.tail(robot.dimv()) = s_new_prev.v - s_old_prev.v;
  dx_ = kkt_matrix_inverse_.forwardCorrectionSerialCoeff() * x_res_;
  robot.integrateConfiguration(dx_.head(robot.dimv()), -1, s_new.q);
  s_new.v.noalias() -= dx_.tail(robot.dimv());
}


void SplitParNMPC::forwardCorrectionParallel(const Robot& robot,
                                             SplitDirection& d,
                                             SplitSolution& s_new) {
  const int dim_kkt = d.dimKKT();
  const int dimx = 2*robot.dimv();
  d.split_direction().segment(0, dim_kkt-dimx) = kkt_matrix_inverse_.forwardCorrectionParallelCoeff() * x_res_;
  s_new.lmd.noalias() -= d.dlmd();
  s_new.gmm.noalias() -= d.dgmm();
  s_new.mu_active().noalias() -= d.dmu();
  s_new.a.noalias() -= d.da();
  s_new.f_active().noalias() -= d.df();
}


void SplitParNMPC::computePrimalAndDualDirection(const Robot& robot, 
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
  id_condenser_.computeCondensedDirection(s, dtau, d);
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, dtau, 
                                             d.da(), d.df(), d.dq(), d.dv(), d.du);
}

 
double SplitParNMPC::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


double SplitParNMPC::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


std::pair<double, double> SplitParNMPC::costAndConstraintsViolation(
    Robot& robot, const double t, const double dtau, const SplitSolution& s,
    const bool is_terminal) {
  assert(dtau > 0);
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, s.q, s.v, s.a, s.f, s.u);
  if (is_terminal) {
    cost += cost_->phi(robot, cost_data_, t, s.q, s.v);
  }
  cost += constraints_->costSlackBarrier(constraints_data_);
  if (robot.has_floating_base()) {
    id_condenser_.getFloatingBaseConstraint(dtau, s, kkt_residual_);
  }
  double constraints_violation = 0;
  constraints_violation += kkt_residual_.Fq().lpNorm<1>();
  constraints_violation += kkt_residual_.Fv().lpNorm<1>();
  constraints_violation += id_condenser_.inverseDynamicsResidualL1Norm(dtau);
  constraints_violation += kkt_residual_.C().lpNorm<1>();
  constraints_violation 
      += constraints_->residualL1Nrom(robot, constraints_data_, dtau, 
                                      s.a, s.f, s.q, s.v, s.u);
  return std::make_pair(cost, constraints_violation);
}


std::pair<double, double> SplitParNMPC::costAndConstraintsViolation(
    Robot& robot, const double step_size, const double t, const double dtau, 
    const Eigen::VectorXd& q_prev, const Eigen::VectorXd& v_prev, 
    const SplitSolution& s, const SplitDirection& d, SplitSolution& s_tmp) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp.q);
  s_tmp.v = s.v + step_size * d.dv();
  s_tmp.a = s.a + step_size * d.da();
  s_tmp.u = s.u + step_size * d.du;
  if (robot.dimf() > 0) {
    s_tmp.f_active() = s.f_active() + step_size * d.df();
    robot.setContactForces(s_tmp.f);
  }
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, s_tmp.q, s_tmp.v,  
                   s_tmp.a, s_tmp.f, s_tmp.u);
  cost += constraints_->costSlackBarrier(constraints_data_, step_size);
  robot.subtractConfiguration(q_prev, s_tmp.q, kkt_residual_.Fq());
  kkt_residual_.Fq().noalias() += dtau * s_tmp.v;
  kkt_residual_.Fv() = v_prev - s_tmp.v + dtau * s_tmp.a;
  if (robot.dimf() > 0) {
    robot.updateKinematics(s_tmp.q, s_tmp.v, s_tmp.a);
    robot.computeBaumgarteResidual(dtau, kkt_residual_.C());
  }
  if (robot.has_floating_base()) {
    id_condenser_.getFloatingBaseConstraint(dtau, s, kkt_residual_);
  }
  double constraints_violation = 0;
  constraints_violation += kkt_residual_.Fq().lpNorm<1>();
  constraints_violation += kkt_residual_.Fv().lpNorm<1>();
  constraints_violation += id_condenser_.inverseDynamicsResidualL1Norm(robot, dtau, s_tmp);
  constraints_violation += kkt_residual_.C().lpNorm<1>();
  constraints_violation 
      += constraints_->residualL1Nrom(robot, constraints_data_, dtau, s_tmp.a, 
                                      s_tmp.f, s_tmp.q, s_tmp.v, s_tmp.u);
  return std::make_pair(cost, constraints_violation);
}


std::pair<double, double> SplitParNMPC::costAndConstraintsViolation(
    Robot& robot, const double step_size, const double t, const double dtau, 
    const Eigen::VectorXd& q_prev, const Eigen::VectorXd& v_prev, 
    const Eigen::VectorXd& dq_prev, const Eigen::VectorXd& dv_prev, 
    const SplitSolution& s, const SplitDirection& d, SplitSolution& s_tmp, 
    const bool is_terminal) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  assert(dq_prev.size() == robot.dimv());
  assert(dv_prev.size() == robot.dimv());
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp.q);
  s_tmp.v = s.v + step_size * d.dv();
  s_tmp.a = s.a + step_size * d.da();
  s_tmp.u = s.u + step_size * d.du;
  if (robot.dimf() > 0) {
    s_tmp.f_active() = s.f_active() + step_size * d.df();
    robot.setContactForces(s_tmp.f);
  }
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, s_tmp.q, s_tmp.v,  
                   s_tmp.a, s_tmp.f, s_tmp.u);
  if (is_terminal) {
    cost += cost_->phi(robot, cost_data_, t, s_tmp.q, s_tmp.v);
  }
  cost += constraints_->costSlackBarrier(constraints_data_, step_size);
  robot.subtractConfiguration(q_prev, s_tmp.q, kkt_residual_.Fq());
  kkt_residual_.Fq().noalias() += dtau * s_tmp.v + step_size * dq_prev;
  kkt_residual_.Fv() = v_prev - s_tmp.v + dtau * s_tmp.a + step_size * dv_prev;
  if (robot.dimf() > 0) {
    robot.updateKinematics(s_tmp.q, s_tmp.v, s_tmp.a);
    robot.computeBaumgarteResidual(dtau, kkt_residual_.C());
  }
  if (robot.has_floating_base()) {
    id_condenser_.getFloatingBaseConstraint(dtau, s, kkt_residual_);
  }
  double constraints_violation = 0;
  constraints_violation += kkt_residual_.Fq().lpNorm<1>();
  constraints_violation += kkt_residual_.Fv().lpNorm<1>();
  constraints_violation += id_condenser_.inverseDynamicsConstraintViolation(robot, dtau, s_tmp);
  constraints_violation += kkt_residual_.C().lpNorm<1>();
  constraints_violation 
      += constraints_->residualL1Nrom(robot, constraints_data_, dtau, s_tmp.a, 
                                      s_tmp.f, s_tmp.q, s_tmp.v, s_tmp.u);
  return std::make_pair(cost, constraints_violation);
}


void SplitParNMPC::updateDual(const double step_size) {
  assert(step_size > 0);
  assert(step_size <= 1);
  constraints_->updateDual(constraints_data_, step_size);
}


void SplitParNMPC::updatePrimal(Robot& robot, const double step_size, 
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


void SplitParNMPC::getStateFeedbackGain(Eigen::MatrixXd& Kq, 
                                        Eigen::MatrixXd& Kv) const {
  // Kq = du_dq_ + du_da_ * Kaq_ + du_df_.leftCols(dimf_) * Kfq_.topRows(dimf_);
  // Kv = du_dv_ + du_da_ * Kav_ + du_df_.leftCols(dimf_) * Kfv_.topRows(dimf_);
}


double SplitParNMPC::squaredKKTErrorNorm(Robot& robot, const double t, 
                                         const double dtau, 
                                         const Eigen::VectorXd& q_prev, 
                                         const Eigen::VectorXd& v_prev, 
                                         const SplitSolution& s,
                                         const Eigen::VectorXd& lmd_next,
                                         const Eigen::VectorXd& gmm_next,
                                         const Eigen::VectorXd& q_next) {
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  assert(lmd_next.size() == robot.dimv());
  assert(gmm_next.size() == robot.dimv());
  assert(q_next.size() == robot.dimq());
  // Clear the KKT matrix and KKT residual.
  kkt_matrix_.setZero();
  kkt_matrix_.setContactStatus(robot);
  kkt_residual_.setZero();
  kkt_residual_.setContactStatus(robot);
  if (robot.dimf() > 0) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  // Forms the Newton residual by linearizing cost, dynamics, and constraints. 
  linearizer_.linearizeCostAndConstraints(robot, cost_, cost_data_, 
                                          constraints_, constraints_data_, 
                                          t, dtau, s, kkt_residual_);
  linearizer_.linearizeStateEquation(robot, dtau, q_prev, v_prev, s, lmd_next, 
                                     gmm_next, q_next, kkt_residual_, kkt_matrix_);
  linearizer_.linearizeContactConstraints(robot, dtau, kkt_residual_, kkt_matrix_);
  // Newton residual with respect to inverse dynamics.
  id_condenser_.setContactStatus(robot);
  id_condenser_.linearizeCostAndConstraints(robot, cost_, cost_data_, 
                                            constraints_, constraints_data_, 
                                            t, dtau, s);
  id_condenser_.linearizeInverseDynamics(robot, s);
  id_condenser_.augmentInverseDynamicsDerivatives(dtau, s, kkt_residual_);
  id_condenser_.getFloatingBaseConstraint(dtau, s, kkt_residual_);
  // Augment the equality constraints 
  kkt_residual_.lq().noalias() += kkt_matrix_.Cq().transpose() * s.mu_active();
  kkt_residual_.lv().noalias() += kkt_matrix_.Cv().transpose() * s.mu_active();
  kkt_residual_.la().noalias() += kkt_matrix_.Ca().transpose() * s.mu_active();
  if (robot.dimf() > 0) {
    kkt_residual_.lf().noalias() += kkt_matrix_.Cf().transpose() * s.mu_active();
  }
  double error = kkt_residual_.squaredKKTErrorNorm();
  error += id_condenser_.squaredKKTErrorNorm(dtau);
  error += constraints_->residualSquaredNrom(robot, constraints_data_, dtau, 
                                             s.a, s.f, s.q, s.v, s.u);
  return error;
}


double SplitParNMPC::squaredKKTErrorNorm(Robot& robot, const double t, 
                                         const double dtau, 
                                         const Eigen::VectorXd& q_prev, 
                                         const Eigen::VectorXd& v_prev, 
                                         const SplitSolution& s) {
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  // Clear the KKT matrix and KKT residual.
  kkt_matrix_.setZero();
  kkt_matrix_.setContactStatus(robot);
  kkt_residual_.setZero();
  kkt_residual_.setContactStatus(robot);
  if (robot.dimf() > 0) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  // Forms the Newton residual by linearizing cost, dynamics, and constraints. 
  linearizer_.linearizeCostAndConstraints(robot, cost_, cost_data_, 
                                          constraints_, constraints_data_, 
                                          t, dtau, s, kkt_residual_);
  linearizer_.linearizeStateEquation(robot, dtau, q_prev, v_prev, s, 
                                     kkt_residual_, kkt_matrix_);
  linearizer_.linearizeContactConstraints(robot, dtau, kkt_residual_, kkt_matrix_);
  // Condense inverse dynamics.
  id_condenser_.setContactStatus(robot);
  id_condenser_.linearizeCostAndConstraints(robot, cost_, cost_data_, 
                                            constraints_, constraints_data_, 
                                            t, dtau, s);
  id_condenser_.linearizeInverseDynamics(robot, s);
  id_condenser_.augmentInverseDynamicsDerivatives(dtau, s, kkt_residual_);
  id_condenser_.getFloatingBaseConstraint(dtau, s, kkt_residual_);
  // Augment the equality constraints 
  kkt_residual_.lq().noalias() += kkt_matrix_.Cq().transpose() * s.mu_active();
  kkt_residual_.lv().noalias() += kkt_matrix_.Cv().transpose() * s.mu_active();
  kkt_residual_.la().noalias() += kkt_matrix_.Ca().transpose() * s.mu_active();
  if (robot.dimf() > 0) {
    kkt_residual_.lf().noalias() += kkt_matrix_.Cf().transpose() * s.mu_active();
  }
  double error = kkt_residual_.squaredKKTErrorNorm();
  error += id_condenser_.squaredKKTErrorNorm(dtau);
  error += constraints_->residualSquaredNrom(robot, constraints_data_, dtau, 
                                             s.a, s.f, s.q, s.v, s.u);
  return error;
}

} // namespace idocp