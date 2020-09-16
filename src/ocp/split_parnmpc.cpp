#include "idocp/ocp/split_parnmpc.hpp"

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
    state_equation_(robot),
    robot_dynamics_(robot),
    dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    dimKKT_(kkt_residual_.dimKKT()),
    kkt_matrix_inverse_(Eigen::MatrixXd::Zero(kkt_residual_.max_dimKKT(), 
                                              kkt_residual_.max_dimKKT())),
    x_res_(Eigen::VectorXd::Zero(2*robot.dimv())),
    dx_(Eigen::VectorXd::Zero(2*robot.dimv())),
    use_kinematics_(false) {
  if (cost_->useKinematics() || constraints_->useKinematics() 
                             || robot.max_dimf() > 0) {
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
    state_equation_(),
    robot_dynamics_(),
    dimv_(0),
    dimx_(0),
    dimKKT_(0),
    kkt_matrix_inverse_(),
    x_res_(),
    dx_(),
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
  constraints_->setSlackAndDual(robot, constraints_data_, dtau, s);
}


void SplitParNMPC::coarseUpdate(Robot& robot, const double t, const double dtau, 
                                const Eigen::VectorXd& q_prev, 
                                const Eigen::VectorXd& v_prev,
                                const SplitSolution& s,
                                const SplitSolution& s_next,
                                const Eigen::MatrixXd& aux_mat_next_old,
                                SplitDirection& d, 
                                SplitSolution& s_new_coarse) {
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  assert(s.dimf() == robot.dimf());
  assert(s.dimc() == robot.dim_passive()+robot.dimf());
  assert(aux_mat_next_old.rows() == 2*robot.dimv());
  assert(aux_mat_next_old.cols() == 2*robot.dimv());
  assert(d.dimf() == robot.dimf());
  assert(d.dimc() == robot.dim_passive()+robot.dimf());
  assert(s_new_coarse.dimf() == robot.dimf());
  assert(s_new_coarse.dimc() == robot.dim_passive()+robot.dimf());
  kkt_residual_.setContactStatus(robot);
  kkt_matrix_.setContactStatus(robot);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  // condensing the inverse dynamics
  kkt_residual_.lu.setZero();
  kkt_matrix_.Quu.setZero();
  cost_->lu(robot, cost_data_, t, dtau, s.u, kkt_residual_.lu);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, 
                                    kkt_residual_.lu);
  cost_->luu(robot, cost_data_, t, dtau, s.u, kkt_matrix_.Quu);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dtau, s.u, 
                                     kkt_matrix_.Quu, kkt_residual_.lu);
  robot_dynamics_.condenseRobotDynamics(robot, dtau, s, kkt_matrix_, 
                                        kkt_residual_);
  // forms the KKT matrix and KKT residual
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, 
                                    kkt_residual_);
  state_equation_.linearizeBackwardEuler(robot, dtau, q_prev, v_prev, s, 
                                         s_next.lmd, s_next.gmm, s_next.q, 
                                         kkt_matrix_, kkt_residual_);
  cost_->computeStageCostHessian(robot, cost_data_, t, dtau, s, kkt_matrix_);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dtau, s, 
                                     kkt_matrix_, kkt_residual_);
  // coarse update of the solution
  kkt_matrix_.Qxx().noalias() += aux_mat_next_old;
  kkt_matrix_.symmetrize(); 
  dimKKT_ = kkt_residual_.dimKKT();
  kkt_matrix_.invert(dtau, kkt_matrix_inverse_.topLeftCorner(dimKKT_, dimKKT_));
  d.split_direction() = kkt_matrix_inverse_.topLeftCorner(dimKKT_, dimKKT_)
                          * kkt_residual_.KKT_residual();
  s_new_coarse.lmd = s.lmd - d.dlmd();
  s_new_coarse.gmm = s.gmm - d.dgmm();
  s_new_coarse.mu_active() = s.mu_active() - d.dmu();
  s_new_coarse.a = s.a - d.da();
  s_new_coarse.f_active() = s.f_active() - d.df();
  robot.integrateConfiguration(s.q, d.dq(), -1, s_new_coarse.q);
  s_new_coarse.v = s.v - d.dv();
}


void SplitParNMPC::coarseUpdateTerminal(Robot& robot, const double t, 
                                        const double dtau, 
                                        const Eigen::VectorXd& q_prev, 
                                        const Eigen::VectorXd& v_prev, 
                                        const SplitSolution& s, 
                                        SplitDirection& d, 
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
  kkt_residual_.setContactStatus(robot);
  kkt_matrix_.setContactStatus(robot);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  // condensing the inverse dynamics
  kkt_residual_.lu.setZero();
  kkt_matrix_.Quu.setZero();
  cost_->lu(robot, cost_data_, t, dtau, s.u, kkt_residual_.lu);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, 
                                    kkt_residual_.lu);
  cost_->luu(robot, cost_data_, t, dtau, s.u, kkt_matrix_.Quu);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dtau, s.u, 
                                     kkt_matrix_.Quu, kkt_residual_.lu);
  robot_dynamics_.condenseRobotDynamics(robot, dtau, s, kkt_matrix_, 
                                        kkt_residual_);
  // forms the KKT matrix and KKT residual
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  cost_->computeTerminalCostDerivatives(robot, cost_data_, t, s, kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, 
                                    kkt_residual_);
  state_equation_.linearizeBackwardEulerTerminal(robot, dtau, q_prev, v_prev, s, 
                                                 kkt_matrix_, kkt_residual_);
  cost_->computeStageCostHessian(robot, cost_data_, t, dtau, s, kkt_matrix_);
  cost_->computeTerminalCostHessian(robot, cost_data_, t, s, kkt_matrix_);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dtau, s, 
                                     kkt_matrix_, kkt_residual_);
  // coarse update of the solution
  dimKKT_ = kkt_residual_.dimKKT();
  kkt_matrix_.symmetrize(); 
  dimKKT_ = kkt_residual_.dimKKT();
  kkt_matrix_.invert(dtau, kkt_matrix_inverse_.topLeftCorner(dimKKT_, dimKKT_));
  d.split_direction() = kkt_matrix_inverse_.topLeftCorner(dimKKT_, dimKKT_)
                          * kkt_residual_.KKT_residual();
  s_new_coarse.lmd = s.lmd - d.dlmd();
  s_new_coarse.gmm = s.gmm - d.dgmm();
  s_new_coarse.mu_active() = s.mu_active() - d.dmu();
  s_new_coarse.a = s.a - d.da();
  s_new_coarse.f_active() = s.f_active() - d.df();
  robot.integrateConfiguration(s.q, d.dq(), -1, s_new_coarse.q);
  s_new_coarse.v = s.v - d.dv();
}


void SplitParNMPC::getAuxiliaryMatrix(Eigen::MatrixXd& aux_mat) const {
  assert(aux_mat.rows() == dimx_);
  assert(aux_mat.cols() == dimx_);
  aux_mat = - kkt_matrix_inverse_.topLeftCorner(dimx_, dimx_);
}


void SplitParNMPC::getTerminalCostHessian(Robot& robot, const double t, 
                                          const SplitSolution& s, 
                                          Eigen::MatrixXd& phixx) {
  assert(phixx.rows() == dimx_);
  assert(phixx.cols() == dimx_);
  kkt_matrix_.Qxx().setZero();
  cost_->computeTerminalCostHessian(robot, cost_data_, t, s, kkt_matrix_);
  phixx = kkt_matrix_.Qxx();
  kkt_matrix_.Qxx().setZero();
}


void SplitParNMPC::backwardCorrectionSerial(const Robot& robot,
                                            const SplitSolution& s_next, 
                                            const SplitSolution& s_new_next,
                                            SplitSolution& s_new) {
  x_res_.head(robot.dimv()) = s_new_next.lmd - s_next.lmd;
  x_res_.tail(robot.dimv()) = s_new_next.gmm - s_next.gmm;
  dx_.noalias() 
      = kkt_matrix_inverse_.block(0, dimKKT_-dimx_, dimx_, dimx_) * x_res_;
  s_new.lmd.noalias() -= dx_.head(robot.dimv());
  s_new.gmm.noalias() -= dx_.tail(robot.dimv());
}


void SplitParNMPC::backwardCorrectionParallel(const Robot& robot, 
                                              SplitDirection& d,
                                              SplitSolution& s_new) {
  d.split_direction().tail(dimKKT_-dimx_).noalias()
      = kkt_matrix_inverse_.block(dimx_, dimKKT_-dimx_, dimKKT_-dimx_, dimx_) 
          * x_res_;
  s_new.mu_active().noalias() -= d.dmu();
  s_new.a.noalias() -= d.da();
  s_new.f_active().noalias() -= d.df();
  robot.integrateConfiguration(d.dq(), -1, s_new.q);
  s_new.v.noalias() -= d.dv();
}


void SplitParNMPC::forwardCorrectionSerial(const Robot& robot,
                                           const SplitSolution& s_prev,
                                           const SplitSolution& s_new_prev, 
                                           SplitSolution& s_new) {
  robot.subtractConfiguration(s_new_prev.q, s_prev.q, 
                              x_res_.head(robot.dimv()));
  x_res_.tail(robot.dimv()) = s_new_prev.v - s_prev.v;
  dx_.noalias() = kkt_matrix_inverse_.block(dimKKT_-dimx_, 0, dimx_, dimx_) 
                    * x_res_;
  robot.integrateConfiguration(dx_.head(robot.dimv()), -1, s_new.q);
  s_new.v.noalias() -= dx_.tail(robot.dimv());
}


void SplitParNMPC::forwardCorrectionParallel(const Robot& robot,
                                             SplitDirection& d,
                                             SplitSolution& s_new) {
  d.split_direction().head(dimKKT_-dimx_).noalias()
      = kkt_matrix_inverse_.topLeftCorner(dimKKT_-dimx_, dimx_) * x_res_;
  s_new.lmd.noalias() -= d.dlmd();
  s_new.gmm.noalias() -= d.dgmm();
  s_new.mu_active().noalias() -= d.dmu();
  s_new.a.noalias() -= d.da();
  s_new.f_active().noalias() -= d.df();
}


void SplitParNMPC::computePrimalAndDualDirection(Robot& robot, 
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
  robot_dynamics_.computeCondensedDirection(dtau, kkt_matrix_, kkt_residual_, d);
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, dtau, d);
}

 
double SplitParNMPC::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


double SplitParNMPC::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


std::pair<double, double> SplitParNMPC::costAndViolation(
    Robot& robot, const double t, const double dtau, const SplitSolution& s) {
  assert(dtau > 0);
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, s);
  cost += constraints_->costSlackBarrier(constraints_data_);
  double violation = 0;
  violation += state_equation_.violationL1Norm(kkt_residual_);
  violation += robot_dynamics_.violationL1Norm(robot, dtau, s, kkt_residual_);
  violation += constraints_->residualL1Nrom(robot, constraints_data_, dtau, s);
  return std::make_pair(cost, violation);
}


std::pair<double, double> SplitParNMPC::costAndViolation(
    Robot& robot, const double step_size, const double t, const double dtau, 
    const Eigen::VectorXd& q_prev, const Eigen::VectorXd& v_prev, 
    const SplitSolution& s, const SplitDirection& d, SplitSolution& s_tmp) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  s_tmp.a = s.a + step_size * d.da();
  if (robot.has_active_contacts()) {
    s_tmp.f_active() = s.f_active() + step_size * d.df();
    robot.setContactForces(s_tmp.f);
  }
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp.q);
  s_tmp.v = s.v + step_size * d.dv();
  s_tmp.u = s.u + step_size * d.du;
  if (use_kinematics_) {
    robot.updateKinematics(s_tmp.q, s_tmp.v, s_tmp.a);
  }
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, s_tmp);
  cost += constraints_->costSlackBarrier(constraints_data_, step_size);
  double violation = 0;
  violation += state_equation_.computeBackwardEulerViolationL1Norm(
      robot, dtau, q_prev, v_prev, s_tmp, kkt_residual_);
  violation += robot_dynamics_.computeViolationL1Norm(robot, dtau, s_tmp, 
                                                      kkt_residual_);
  violation += constraints_->residualL1Nrom(robot, constraints_data_, dtau, 
                                            s_tmp);
  return std::make_pair(cost, violation);
}


std::pair<double, double> SplitParNMPC::costAndViolation(
    Robot& robot, const double step_size, const double t, const double dtau, 
    const SplitSolution& s_prev, const SplitDirection& d_prev, 
    const SplitSolution& s, const SplitDirection& d, SplitSolution& s_tmp) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  s_tmp.a = s.a + step_size * d.da();
  if (robot.has_active_contacts()) {
    s_tmp.f_active() = s.f_active() + step_size * d.df();
    robot.setContactForces(s_tmp.f);
  }
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp.q);
  s_tmp.v = s.v + step_size * d.dv();
  s_tmp.u = s.u + step_size * d.du;
  if (use_kinematics_) {
    robot.updateKinematics(s_tmp.q, s_tmp.v, s_tmp.a);
  }
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, s_tmp);
  cost += constraints_->costSlackBarrier(constraints_data_, step_size);
  double violation = 0;
  violation += state_equation_.computeBackwardEulerViolationL1Norm(
      robot, step_size, dtau, s_prev.q, s_prev.v, d_prev.dq(), d_prev.dv(), 
      s_tmp, kkt_residual_);
  violation += robot_dynamics_.computeViolationL1Norm(robot, dtau, s_tmp, 
                                                      kkt_residual_);
  violation += constraints_->residualL1Nrom(robot, constraints_data_, dtau, 
                                            s_tmp);
  return std::make_pair(cost, violation);
}


std::pair<double, double> SplitParNMPC::costAndViolationTerminal(
    Robot& robot, const double t, const double dtau, const SplitSolution& s) {
  assert(dtau > 0);
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, s);
  cost += cost_->phi(robot, cost_data_, t, s);
  cost += constraints_->costSlackBarrier(constraints_data_);
  double violation = 0;
  violation += state_equation_.violationL1Norm(kkt_residual_);
  violation += robot_dynamics_.violationL1Norm(robot, dtau, s, kkt_residual_);
  violation += constraints_->residualL1Nrom(robot, constraints_data_, dtau, s);
  return std::make_pair(cost, violation);
}


std::pair<double, double> SplitParNMPC::costAndViolationTerminal(
    Robot& robot, const double step_size, const double t, const double dtau, 
    const SplitSolution& s_prev, const SplitDirection& d_prev,
    const SplitSolution& s, const SplitDirection& d, SplitSolution& s_tmp) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  s_tmp.a = s.a + step_size * d.da();
  if (robot.has_active_contacts()) {
    s_tmp.f_active() = s.f_active() + step_size * d.df();
    robot.setContactForces(s_tmp.f);
  }
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp.q);
  s_tmp.v = s.v + step_size * d.dv();
  s_tmp.u = s.u + step_size * d.du;
  if (use_kinematics_) {
    robot.updateKinematics(s_tmp.q, s_tmp.v, s_tmp.a);
  }
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, s_tmp);
  cost += cost_->phi(robot, cost_data_, t, s_tmp);
  cost += constraints_->costSlackBarrier(constraints_data_, step_size);
  double violation = 0;
  violation += state_equation_.computeBackwardEulerViolationL1Norm(
      robot, step_size, dtau, s_prev.q, s_prev.v, d_prev.dq(), d_prev.dv(), 
      s_tmp, kkt_residual_);
  violation += robot_dynamics_.computeViolationL1Norm(robot, dtau, s_tmp, 
                                                      kkt_residual_);
  violation += constraints_->residualL1Nrom(robot, constraints_data_, dtau, 
                                            s_tmp);
  return std::make_pair(cost, violation);
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
  assert(Kq.rows() == dimv_);
  assert(Kq.cols() == dimv_);
  assert(Kv.rows() == dimv_);
  assert(Kv.cols() == dimv_);
  const int dimc = kkt_residual_.dimc();
  const int dimf = kkt_residual_.dimf();
  const int a_begin = dimx_ + dimc;
  const int f_begin = a_begin + dimv_;
  const int q_begin = f_begin + dimf;
  const int v_begin = q_begin + dimv_;
  robot_dynamics_.getControlInputTorquesSensitivitiesWithRespectToState(
      kkt_matrix_inverse_.block(a_begin, q_begin, dimv_, dimv_),
      kkt_matrix_inverse_.block(a_begin, v_begin, dimv_, dimv_),
      kkt_matrix_inverse_.block(f_begin, q_begin, dimf, dimv_),
      kkt_matrix_inverse_.block(f_begin, v_begin, dimf, dimv_), Kq, Kv);
}


double SplitParNMPC::condensedSquaredKKTErrorNorm(Robot& robot, const double t, 
                                                  const double dtau, 
                                                  const SplitSolution& s) {
  assert(dtau > 0);
  double error = kkt_residual_.KKT_residual().squaredNorm();
  error += constraints_->squaredKKTErrorNorm(robot, constraints_data_, dtau, s);
  return error;
}


double SplitParNMPC::condensedSquaredKKTErrorNormTerminal(
    Robot& robot, const double t, const double dtau, const SplitSolution& s) {
  assert(dtau > 0);
  return condensedSquaredKKTErrorNorm(robot, t, dtau, s);
}


double SplitParNMPC::computeSquaredKKTErrorNorm(Robot& robot, const double t, 
                                                const double dtau, 
                                                const Eigen::VectorXd& q_prev, 
                                                const Eigen::VectorXd& v_prev, 
                                                const SplitSolution& s,
                                                const SplitSolution& s_next) {
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  kkt_matrix_.setContactStatus(robot);
  kkt_residual_.setContactStatus(robot);
  kkt_residual_.setZeroMinimum();
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  cost_->lu(robot, cost_data_, t, dtau, s.u, kkt_residual_.lu);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, 
                                    kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, 
                                    kkt_residual_.lu);
  state_equation_.linearizeBackwardEuler(robot, dtau, q_prev, v_prev, s, 
                                         s_next.lmd, s_next.gmm, s_next.q, 
                                         kkt_matrix_, kkt_residual_);
  robot_dynamics_.augmentRobotDynamics(robot, dtau, s, kkt_matrix_, 
                                       kkt_residual_);
  double error = kkt_residual_.squaredKKTErrorNorm(dtau);
  error += constraints_->squaredKKTErrorNorm(robot, constraints_data_, dtau, s);
  return error;
}


double SplitParNMPC::computeSquaredKKTErrorNormTerminal(
    Robot& robot, const double t, const double dtau, 
    const Eigen::VectorXd& q_prev, const Eigen::VectorXd& v_prev, 
    const SplitSolution& s) {
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  kkt_matrix_.setContactStatus(robot);
  kkt_residual_.setContactStatus(robot);
  kkt_residual_.setZeroMinimum();
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  cost_->lu(robot, cost_data_, t, dtau, s.u, kkt_residual_.lu);
  cost_->computeTerminalCostDerivatives(robot, cost_data_, t, s, kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, 
                                    kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, 
                                    kkt_residual_.lu);
  state_equation_.linearizeBackwardEulerTerminal(robot, dtau, q_prev, v_prev, s, 
                                                 kkt_matrix_, kkt_residual_);
  robot_dynamics_.augmentRobotDynamics(robot, dtau, s, kkt_matrix_, 
                                       kkt_residual_);
  double error = kkt_residual_.squaredKKTErrorNorm(dtau);
  error += constraints_->squaredKKTErrorNorm(robot, constraints_data_, dtau, s);
  return error;
}

} // namespace idocp