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


void SplitParNMPC::coarseUpdate(Robot& robot, 
                                const ContactStatus& contact_status,
                                const double t, const double dtau, 
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
  assert(aux_mat_next_old.rows() == 2*robot.dimv());
  assert(aux_mat_next_old.cols() == 2*robot.dimv());
  setContactStatusForKKT(contact_status);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  // condensing the inverse dynamics
  kkt_residual_.lu.setZero();
  kkt_matrix_.Quu.setZero();
  cost_->lu(robot, cost_data_, t, dtau, s.u, kkt_residual_.lu);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s.u,
                                    kkt_residual_.lu);
  cost_->luu(robot, cost_data_, t, dtau, s.u, kkt_matrix_.Quu);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dtau, s.u, 
                                     kkt_matrix_.Quu, kkt_residual_.lu);
  robot_dynamics_.condenseRobotDynamics(robot, contact_status, dtau, s, 
                                        kkt_matrix_, kkt_residual_);
  // construct the KKT matrix and KKT residual
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s,
                                    kkt_residual_);
  stateequation::LinearizeBackwardEuler(robot, dtau, q_prev, v_prev, s, s_next, 
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
  s_new_coarse.mu_stack() = s.mu_stack() - d.dmu();
  s_new_coarse.a = s.a - d.da();
  s_new_coarse.f_stack() = s.f_stack() - d.df();
  robot.integrateConfiguration(s.q, d.dq(), -1, s_new_coarse.q);
  s_new_coarse.v = s.v - d.dv();
}


void SplitParNMPC::coarseUpdateTerminal(Robot& robot,  
                                        const ContactStatus& contact_status,
                                        const double t, const double dtau, 
                                        const Eigen::VectorXd& q_prev, 
                                        const Eigen::VectorXd& v_prev, 
                                        const SplitSolution& s, 
                                        SplitDirection& d, 
                                        SplitSolution& s_new_coarse) {
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  setContactStatusForKKT(contact_status);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  // condensing the inverse dynamics
  kkt_residual_.lu.setZero();
  kkt_matrix_.Quu.setZero();
  cost_->lu(robot, cost_data_, t, dtau, s.u, kkt_residual_.lu);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s.u,
                                    kkt_residual_.lu);
  cost_->luu(robot, cost_data_, t, dtau, s.u, kkt_matrix_.Quu);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dtau, s.u, 
                                     kkt_matrix_.Quu, kkt_residual_.lu);
  robot_dynamics_.condenseRobotDynamics(robot, contact_status, dtau, s, 
                                        kkt_matrix_, kkt_residual_);
  // constructs the KKT matrix and KKT residual
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  cost_->computeTerminalCostDerivatives(robot, cost_data_, t, s, kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s,
                                    kkt_residual_);
  stateequation::LinearizeBackwardEulerTerminal(robot, dtau, q_prev, v_prev, s, 
                                                kkt_matrix_, kkt_residual_);
  cost_->computeStageCostHessian(robot, cost_data_, t, dtau, s, kkt_matrix_);
  cost_->computeTerminalCostHessian(robot, cost_data_, t, s, kkt_matrix_);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dtau, s, 
                                     kkt_matrix_, kkt_residual_);
  // coarse update of the solution
  kkt_matrix_.symmetrize(); 
  dimKKT_ = kkt_residual_.dimKKT();
  kkt_matrix_.invert(dtau, kkt_matrix_inverse_.topLeftCorner(dimKKT_, dimKKT_));
  d.split_direction() = kkt_matrix_inverse_.topLeftCorner(dimKKT_, dimKKT_)
                          * kkt_residual_.KKT_residual();
  s_new_coarse.lmd = s.lmd - d.dlmd();
  s_new_coarse.gmm = s.gmm - d.dgmm();
  s_new_coarse.mu_stack() = s.mu_stack() - d.dmu();
  s_new_coarse.a = s.a - d.da();
  s_new_coarse.f_stack() = s.f_stack() - d.df();
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
                                              SplitSolution& s_new) const {
  d.split_direction().tail(dimKKT_-dimx_).noalias()
      = kkt_matrix_inverse_.block(dimx_, dimKKT_-dimx_, dimKKT_-dimx_, dimx_) 
          * x_res_;
  s_new.mu_stack().noalias() -= d.dmu();
  s_new.a.noalias() -= d.da();
  s_new.f_stack().noalias() -= d.df();
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
                                             SplitSolution& s_new) const {
  d.split_direction().head(dimKKT_-dimx_).noalias()
      = kkt_matrix_inverse_.topLeftCorner(dimKKT_-dimx_, dimx_) * x_res_;
  s_new.lmd.noalias() -= d.dlmd();
  s_new.gmm.noalias() -= d.dgmm();
  s_new.mu_stack().noalias() -= d.dmu();
  s_new.a.noalias() -= d.da();
  s_new.f_stack().noalias() -= d.df();
}


void SplitParNMPC::computePrimalAndDualDirection(Robot& robot, 
                                                 const double dtau, 
                                                 const SplitSolution& s, 
                                                 const SplitSolution& s_new,
                                                 SplitDirection& d) {
  d.dlmd() = s_new.lmd - s.lmd;
  d.dgmm() = s_new.gmm - s.gmm;
  d.dmu() = s_new.mu_stack() - s.mu_stack();
  d.da() = s_new.a - s.a;
  d.df() = s_new.f_stack() - s.f_stack();
  robot.subtractConfiguration(s_new.q, s.q, d.dq());
  d.dv() = s_new.v - s.v;
  robot_dynamics_.computeCondensedDirection(dtau, kkt_matrix_, kkt_residual_, d);
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, dtau, s, d);
}

 
double SplitParNMPC::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


double SplitParNMPC::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


std::pair<double, double> SplitParNMPC::costAndConstraintViolation(
    Robot& robot, const double t, const double dtau, const SplitSolution& s) {
  assert(dtau > 0);
  return std::make_pair(cost(robot, t, dtau, s), constraintViolation(dtau));
}


std::pair<double, double> SplitParNMPC::costAndConstraintViolation(
    Robot& robot, const ContactStatus& contact_status, const double step_size, 
    const double t, const double dtau, const Eigen::VectorXd& q_prev, 
    const Eigen::VectorXd& v_prev, const SplitSolution& s, 
    const SplitDirection& d, SplitSolution& s_tmp) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  setContactStatusForKKT(contact_status);
  s_tmp.setContactStatus(contact_status);
  s_tmp.a = s.a + step_size * d.da();
  if (contact_status.hasActiveContacts()) {
    s_tmp.setContactStatus(contact_status);
    s_tmp.f_stack() = s.f_stack() + step_size * d.df();
    s_tmp.set_f();
    robot.setContactForces(contact_status, s_tmp.f);
  }
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp.q);
  s_tmp.v = s.v + step_size * d.dv();
  s_tmp.u = s.u + step_size * d.du;
  if (use_kinematics_) {
    robot.updateKinematics(s_tmp.q, s_tmp.v, s_tmp.a);
  }
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, dtau, 
                                             s_tmp);
  stateequation::ComputeBackwardEulerResidual(robot, dtau, q_prev, v_prev, 
                                              s_tmp, kkt_residual_);
  robot_dynamics_.computeRobotDynamicsResidual(robot, contact_status, dtau, 
                                               s_tmp, kkt_residual_);
  return std::make_pair(cost(robot, t, dtau, s_tmp), constraintViolation(dtau));
}


std::pair<double, double> SplitParNMPC::costAndConstraintViolation(
    Robot& robot, const ContactStatus& contact_status, const double step_size, 
    const double t, const double dtau, const SplitSolution& s_prev, 
    const SplitDirection& d_prev, const SplitSolution& s, 
    const SplitDirection& d, SplitSolution& s_tmp) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  setContactStatusForKKT(contact_status);
  s_tmp.setContactStatus(contact_status);
  s_tmp.a = s.a + step_size * d.da();
  if (contact_status.hasActiveContacts()) {
    s_tmp.setContactStatus(contact_status);
    s_tmp.f_stack() = s.f_stack() + step_size * d.df();
    s_tmp.set_f();
    robot.setContactForces(contact_status, s_tmp.f);
  }
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp.q);
  s_tmp.v = s.v + step_size * d.dv();
  s_tmp.u = s.u + step_size * d.du;
  if (use_kinematics_) {
    robot.updateKinematics(s_tmp.q, s_tmp.v, s_tmp.a);
  }
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, dtau, 
                                             s_tmp);
  stateequation::ComputeBackwardEulerResidual(robot, step_size, dtau, s_prev.q, 
                                              s_prev.v, d_prev.dq(), d_prev.dv(), 
                                              s_tmp, kkt_residual_);
  robot_dynamics_.computeRobotDynamicsResidual(robot, contact_status, dtau, 
                                               s_tmp, kkt_residual_);
  return std::make_pair(cost(robot, t, dtau, s_tmp), constraintViolation(dtau));
}


std::pair<double, double> SplitParNMPC::costAndConstraintViolationTerminal(
    Robot& robot, const double t, const double dtau, const SplitSolution& s) {
  assert(dtau > 0);
  return std::make_pair(costTerminal(robot, t, dtau, s), 
                        constraintViolation(dtau));
}


std::pair<double, double> SplitParNMPC::costAndConstraintViolationTerminal(
    Robot& robot, const ContactStatus& contact_status, const double step_size, 
    const double t, const double dtau, const SplitSolution& s_prev, 
    const SplitDirection& d_prev, const SplitSolution& s, 
    const SplitDirection& d, SplitSolution& s_tmp) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  setContactStatusForKKT(contact_status);
  s_tmp.setContactStatus(contact_status);
  s_tmp.a = s.a + step_size * d.da();
  if (contact_status.hasActiveContacts()) {
    s_tmp.f_stack() = s.f_stack() + step_size * d.df();
    s_tmp.set_f();
    robot.setContactForces(contact_status, s_tmp.f);
  }
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp.q);
  s_tmp.v = s.v + step_size * d.dv();
  s_tmp.u = s.u + step_size * d.du;
  if (use_kinematics_) {
    robot.updateKinematics(s_tmp.q, s_tmp.v, s_tmp.a);
  }
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, dtau, 
                                             s_tmp);
  stateequation::ComputeBackwardEulerResidual(robot, step_size, dtau, s_prev.q, 
                                              s_prev.v, d_prev.dq(), d_prev.dv(), 
                                              s_tmp, kkt_residual_);
  robot_dynamics_.computeRobotDynamicsResidual(robot, contact_status, dtau, 
                                               s_tmp, kkt_residual_);
  return std::make_pair(costTerminal(robot, t, dtau, s_tmp), 
                        constraintViolation(dtau));
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
  s.mu_stack().noalias() += step_size * d.dmu();
  s.set_mu_contact();
  s.a.noalias() += step_size * d.da();
  s.f_stack().noalias() += step_size * d.df();
  s.set_f();
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
  robot_dynamics_.getStateFeedbackGain(
      kkt_matrix_inverse_.block(a_begin, q_begin, dimv_, dimv_),
      kkt_matrix_inverse_.block(a_begin, v_begin, dimv_, dimv_),
      kkt_matrix_inverse_.block(f_begin, q_begin, dimf, dimv_),
      kkt_matrix_inverse_.block(f_begin, v_begin, dimf, dimv_), Kq, Kv);
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
  kkt_residual_.setZeroMinimum();
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  cost_->lu(robot, cost_data_, t, dtau, s.u, kkt_residual_.lu);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s,
                                    kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s.u,
                                    kkt_residual_.lu);
  stateequation::LinearizeBackwardEuler(robot, dtau, q_prev, v_prev, s, s_next,
                                        kkt_matrix_, kkt_residual_);
  robot_dynamics_.linearizeRobotDynamics(robot, contact_status, dtau, s, 
                                         kkt_matrix_, kkt_residual_);
}


void SplitParNMPC::computeKKTResidualTerminal(
    Robot& robot, const ContactStatus& contact_status, const double t, 
    const double dtau, const Eigen::VectorXd& q_prev, 
    const Eigen::VectorXd& v_prev, const SplitSolution& s) {
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  setContactStatusForKKT(contact_status);
  kkt_residual_.setZeroMinimum();
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  cost_->lu(robot, cost_data_, t, dtau, s.u, kkt_residual_.lu);
  cost_->computeTerminalCostDerivatives(robot, cost_data_, t, s, kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s,
                                    kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s.u,
                                    kkt_residual_.lu);
  stateequation::LinearizeBackwardEulerTerminal(robot, dtau, q_prev, v_prev, s, 
                                                kkt_matrix_, kkt_residual_);
  robot_dynamics_.linearizeRobotDynamics(robot, contact_status, dtau, s, 
                                         kkt_matrix_, kkt_residual_);
}


double SplitParNMPC::squaredNormKKTResidual(const double dtau) const {
  double error = 0;
  error += kkt_residual_.lq().squaredNorm();
  error += kkt_residual_.lv().squaredNorm();
  error += kkt_residual_.la().squaredNorm();
  error += kkt_residual_.lf().squaredNorm();
  error += kkt_residual_.lu.squaredNorm();
  error += stateequation::SquaredNormStateEuqationResidual(kkt_residual_);
  error += robot_dynamics_.squaredNormRobotDynamicsResidual(dtau, kkt_residual_);
  error += constraints_->squaredNormPrimalAndDualResidual(constraints_data_);
  return error;
}


double SplitParNMPC::cost(Robot& robot, const double t, const double dtau, 
                          const SplitSolution& s) {
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, s);
  cost += constraints_->costSlackBarrier(constraints_data_);
  return cost;
}


double SplitParNMPC::costTerminal(Robot& robot, const double t, 
                                  const double dtau, const SplitSolution& s) {
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, s);
  cost += cost_->phi(robot, cost_data_, t, s);
  cost += constraints_->costSlackBarrier(constraints_data_);
  return cost;
}


double SplitParNMPC::constraintViolation(const double dtau) const {
  double violation = 0;
  violation += stateequation::L1NormStateEuqationResidual(kkt_residual_);
  violation += robot_dynamics_.l1NormRobotDynamicsResidual(dtau, kkt_residual_);
  violation += constraints_->l1NormPrimalResidual(constraints_data_);
  return violation;
}

} // namespace idocp