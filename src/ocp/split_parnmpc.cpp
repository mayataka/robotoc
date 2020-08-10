#include "idocp/ocp/split_parnmpc.hpp"
#include "idocp/ocp/parnmpc_linearizer.hpp"

#include <assert.h>


namespace idocp {

SplitParNMPC::SplitParNMPC(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost,
    const std::shared_ptr<ConstraintsInterface>& constraints) 
  : cost_(cost),
    constraints_(constraints),
    cost_data_(robot),
    joint_constraints_(robot),
    kkt_residual_(robot),
    kkt_matrix_(robot),
    kkt_composition_(robot),
    has_floating_base_(robot.has_floating_base()),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dim_passive_(robot.dim_passive()),
    max_dimf_(robot.max_dimf()),
    max_dimc_(robot.dim_passive()+robot.max_dimf()),
    lu_(Eigen::VectorXd::Zero(robot.dimv())),
    lu_condensed_(Eigen::VectorXd::Zero(robot.dimv())),
    u_res_(Eigen::VectorXd::Zero(robot.dimv())),
    du_(Eigen::VectorXd::Zero(robot.dimv())),
    dq_res_(Eigen::VectorXd::Zero(robot.dimv())),
    dv_res_(Eigen::VectorXd::Zero(robot.dimv())),
    kkt_matrix_inverse_(kkt_matrix_.max_dimKKT(), kkt_matrix_.max_dimKKT()),
    luu_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    du_dq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    du_dv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    du_da_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    du_df_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    dsubtract_dqminus_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dsubtract_dqplus_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    // The following variables are only needed for line search
    q_tmp_(Eigen::VectorXd::Zero(robot.dimq())), 
    v_tmp_(Eigen::VectorXd::Zero(robot.dimv())), 
    a_tmp_(Eigen::VectorXd::Zero(robot.dimv())), 
    f_tmp_(Eigen::VectorXd::Zero(robot.max_dimf())), 
    u_tmp_(Eigen::VectorXd::Zero(robot.dimv())), 
    u_res_tmp_(Eigen::VectorXd::Zero(robot.dimv())) {
}


SplitParNMPC::SplitParNMPC() 
  : cost_(),
    constraints_(),
    cost_data_(),
    joint_constraints_(),
    kkt_residual_(),
    kkt_matrix_(),
    kkt_composition_(),
    has_floating_base_(false),
    dimq_(0),
    dimv_(0),
    dim_passive_(0),
    max_dimf_(0),
    max_dimc_(0),
    lu_(),
    lu_condensed_(),
    u_res_(),
    du_(),
    dq_res_(),
    dv_res_(),
    luu_(),
    kkt_matrix_inverse_(),
    du_dq_(),
    du_dv_(),
    du_da_(),
    du_df_(), 
    dsubtract_dqminus_(),
    dsubtract_dqplus_(),
    q_tmp_(), 
    v_tmp_(), 
    a_tmp_(), 
    f_tmp_(), 
    u_tmp_(), 
    u_res_tmp_() {
}


SplitParNMPC::~SplitParNMPC() {
}


bool SplitParNMPC::isFeasible(const Robot& robot, const Eigen::VectorXd& q, 
                              const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                              const Eigen::VectorXd& u) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(u.size() == dimv_);
  return joint_constraints_.isFeasible(q, v, a, u);
}


void SplitParNMPC::initConstraints(const Robot& robot, const int time_step, 
                                   const double dtau, const Eigen::VectorXd& q, 
                                   const Eigen::VectorXd& v, 
                                   const Eigen::VectorXd& a, 
                                   const Eigen::VectorXd& u) {
  assert(time_step >= 0);
  assert(dtau > 0);
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(u.size() == dimv_);
  joint_constraints_.setTimeStep(time_step);
  joint_constraints_.setSlackAndDual(dtau, q, v, a, u);
}


void SplitParNMPC::coarseUpdate(Robot& robot, const double t, const double dtau, 
                                const Eigen::VectorXd& q_prev, 
                                const Eigen::VectorXd& v_prev,
                                const Eigen::VectorXd& lmd, 
                                const Eigen::VectorXd& gmm, 
                                const Eigen::VectorXd& mu, 
                                const Eigen::VectorXd& a,
                                const Eigen::VectorXd& f, 
                                const Eigen::VectorXd& q, 
                                const Eigen::VectorXd& v, 
                                const Eigen::VectorXd& u, 
                                const Eigen::VectorXd& lmd_next,
                                const Eigen::VectorXd& gmm_next,
                                const Eigen::VectorXd& q_next,
                                const Eigen::VectorXd& v_next,
                                const Eigen::MatrixXd& aux_mat_next_old,
                                Eigen::MatrixXd& aux_mat) {
  assert(dtau > 0);
  assert(q_prev.size() == dimq_);
  assert(v_prev.size() == dimv_);
  assert(lmd.size() == dimv_);
  assert(gmm.size() == dimv_);
  assert(mu.size() == max_dimc_);
  assert(a.size() == dimv_);
  assert(f.size() == max_dimf_);
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(u.size() == dimv_);
  assert(lmd_next.size() == dimv_);
  assert(gmm_next.size() == dimv_);
  assert(q_next.size() == dimq_);
  assert(v_next.size() == dimv_);
  assert(aux_mat_next_old.rows() == 2*dimv_);
  assert(aux_mat_next_old.cols() == 2*dimv_);
  assert(aux_mat.rows() == 2*dimv_);
  assert(aux_mat.cols() == 2*dimv_);
  // Reset the KKT matrix and KKT residual.
  kkt_matrix_.setZero();
  kkt_matrix_.setContactStatus(robot);
  kkt_residual_.setZero();
  kkt_residual_.setContactStatus(robot);
  // Extract active blocks for matrices with associated with contacts.
  const Eigen::Ref<const Eigen::MatrixXd> du_df_active 
      = du_df_.leftCols(robot.dimf());
  const Eigen::Ref<const Eigen::VectorXd> mu_active 
      = mu.head(robot.dim_passive()+robot.dimf());
  const Eigen::Ref<const Eigen::VectorXd> mu_passive 
      = mu.head(robot.dim_passive());
  const Eigen::Ref<const Eigen::VectorXd> f_active = f.head(robot.dimf());
  Eigen::Ref<Eigen::VectorXd> mu_tmp_acitve = mu_tmp_.head(robot.dim_passive());
  Eigen::Ref<Eigen::VectorXd> f_tmp_active = f_tmp_.head(robot.dimf());
  // Compute KKT residual and KKT matrix.
  if (robot.dimf() > 0) {
    robot.updateKinematics(q, v, a);
  }
  parnmpclinearizer::linearizeStageCost(robot, cost_, cost_data_, t, dtau, 
                                        q, v, a, f, u, kkt_residual_, lu_);
  parnmpclinearizer::linearizeDynamics(robot, dtau, q, v, a, f, u,  
                                       q_prev, v_prev, kkt_residual_, u_res_, 
                                       du_dq_, du_dv_, du_da_, du_df_);
  parnmpclinearizer::linearizeConstraints(robot, dtau, u, u_res_,  
                                          du_dq_, du_dv_, du_da_, du_df_, 
                                          kkt_residual_, kkt_matrix_);
  // Condense the control input torques and the Lagrange multiplier with 
  // respect to inverse dynamics.
  joint_constraints_.augmentDualResidual(dtau, lu_);
  cost_->luu(robot, cost_data_, t, dtau, u, luu_);
  joint_constraints_.condenseSlackAndDual(dtau, u, luu_, lu_);
  lu_condensed_ = lu_ + luu_ * u_res_;
  // Augment the condensed Newton residual of the contorl input torques. 
  kkt_residual_.lq().noalias() += du_dq_.transpose() * lu_condensed_;
  kkt_residual_.lv().noalias() += du_dv_.transpose() * lu_condensed_;
  kkt_residual_.la().noalias() += du_da_.transpose() * lu_condensed_;
  if (robot.dimf() > 0) {
    kkt_residual_.lf().noalias() += du_df_active.transpose() * lu_condensed_;
  }
  // Augmnet the partial derivatives of the state equation.
  if (robot.has_floating_base()) {
    robot.dSubtractdConfigurationMinus(q_prev, q, dsubtract_dqminus_);
    robot.dSubtractdConfigurationPlus(q, q_next, dsubtract_dqplus_);
    kkt_residual_.lq().noalias() 
        += dsubtract_dqplus_.transpose() * lmd_next
            + dsubtract_dqminus_.transpose() * lmd;
    kkt_residual_.lv().noalias() += dtau * lmd - gmm + gmm_next;
    kkt_residual_.la().noalias() += dtau * gmm;
  }
  else {
    kkt_residual_.lq().noalias() += lmd_next - lmd;
    kkt_residual_.lv().noalias() += dtau * lmd - gmm + gmm_next;
    kkt_residual_.la().noalias() += dtau * gmm;
  }
  // Augmnet the partial derivatives of the inequality constriants.
  joint_constraints_.augmentDualResidual(dtau, kkt_residual_.lq(), 
                                         kkt_residual_.lv(), 
                                         kkt_residual_.la());
  // Augment the equality constraints 
  kkt_residual_.lq().noalias() += kkt_matrix_.Cq().transpose() * mu_active;
  kkt_residual_.lv().noalias() += kkt_matrix_.Cv().transpose() * mu_active;
  kkt_residual_.la().noalias() += kkt_matrix_.Ca().transpose() * mu_active;
  if (dimf_ > 0) {
    kkt_residual_.lf().noalias() += kkt_matrix_.Cf().transpose() * mu_passive;
  }
  // Augment the condensed Hessian of the contorl input torques. 
  kkt_matrix_.Qaa() = du_da_.transpose() * luu_ * du_da_;
  if (robot.dimf() > 0) {
    kkt_matrix_.Qaf() = du_da_.transpose() * luu_ * du_df_active; 
  }
  kkt_matrix_.Qaq() = du_da_.transpose() * luu_ * du_dq_;
  kkt_matrix_.Qav() = du_da_.transpose() * luu_ * du_dv_;
  if (robot.dimf() > 0) {
    kkt_matrix_.Qff() = du_df_active.transpose() * luu_ * du_df_active;
    kkt_matrix_.Qfq() = du_df_active.transpose() * luu_ * du_dq_;
    kkt_matrix_.Qfv() = du_df_active.transpose() * luu_ * du_dv_;
    kkt_matrix_.Qff() = du_df_active.transpose() * luu_ * du_df_active;
  }
  kkt_matrix_.Qqq() = du_dq_.transpose() * luu_ * du_dq_;
  kkt_matrix_.Qqv() = du_dq_.transpose() * luu_ * du_dv_;
  kkt_matrix_.Qvv() = du_dv_.transpose() * luu_ * du_dv_;
  // Condense the slack and dual variables of the inequality constraints on 
  // the configuration, velocity, and acceleration.
  joint_constraints_.condenseSlackAndDual(dtau, q, v, a, kkt_matrix_.Qqq(),  
                                          kkt_matrix_.Qvv(), kkt_matrix_.Qaa(), 
                                          lq_, lv_, la_);
  // Augment the cost function Hessian. 
  cost_->augment_lqq(robot, cost_data_, t, dtau, q, v, a, kkt_matrix_.Qqq());
  cost_->augment_lvv(robot, cost_data_, t, dtau, q, v, a, kkt_matrix_.Qvv());
  cost_->augment_laa(robot, cost_data_, t, dtau, q, v, a, kkt_matrix_.Qaa());
  if (robot.dimf() > 0) {
    cost_->augment_lff(robot, cost_data_, t, dtau, f, kkt_matrix_.Qff());
  }
  kkt_matrix_.symmetrize();
  kkt_matrix_.Qxx().noalias() += aux_mat_next_old;
  const int dim_kkt = kkt_matrix_.dimKKT();
  kkt_matrix_.invert(kkt_matrix_inverse_.topLeftCorner(dim_kkt, dim_kkt));
  aux_mat = - kkt_matrix_inverse_.topLeftCorner(2*robot.dimv(), 2*robot.dimv());
  // coarse update of the solution
  kkt_composition_.set(robot);
  dkkt_.topLeftCorner(dim_kkt, dim_kkt) 
      = kkt_matrix_inverse_.topLeftCorner(dim_kkt, dim_kkt) 
          * kkt_residual.KKT_residual();
  lmd_tmp_ = lmd - dkkt_.segment(kkt_composition_.Fq_begin(), 
                                 kkt_composition_.Fq_size());
  gmm_tmp_ = gmm - dkkt_.segment(kkt_composition_.Fv_begin(), 
                                 kkt_composition_.Fv_size());
  mu_tmp_active = mu_active - dkkt_.segment(kkt_composition_.C_begin(), 
                                            kkt_composition_.C_size());
  a_tmp_ = a - dkkt_.segment(kkt_composition_.Qa_begin(), 
                             kkt_composition_.Qa_size());
  f_tmp_active = f_active - dkkt_.segment(kkt_composition_.Qf_begin(), 
                                          kkt_composition_.Qf_size());
  q_tmp_ = q;
  robot.integrateConfiguration(
      dkkt_.segment(kkt_composition_.Qq_begin(), kkt_composition_.Qq_size()), 
      -1, q_tmp_);
  v_tmp_ = v - dkkt_.segment(kkt_composition_.Qv_begin(), 
                             kkt_composition_.Qv_size());
}


void SplitParNMPC::backwardCollectionSerial(const Eigen::VectorXd& lmd_next_old, 
                                            const Eigen::VectorXd& gmm_next_old, 
                                            const Eigen::VectorXd& lmd_next, 
                                            const Eigen::VectorXd& gmm_next,
                                            Eigen::VectorXd& lmd, 
                                            Eigen::VectorXd& gmm) {
  x_res_.head(dimv_) = lmd_next - lmd_next_old;
  x_res_.tail(dimv_) = gmm_next - gmm_next_old;
  dx_ = kkt_matrix_inverse_.block(kkt_composition_.Fx_begin(), 
                                  kkt_composition_.Qx_begin(), 
                                  kkt_composition_.Fx_size(), 
                                  kkt_composition_.Qx_size()) * x_res_;
  lmd.noalias() -= dx_.head(dimv_);
  gmm.noalias() -= dx_.tail(dimv_);
}


void SplitParNMPC::backwardCollectionParallel(const Robot& robot) {
  const int dim_kkt = kkt_matrix_.dimKKT();
  const int dimx = 2*robot.dimv();
  dkkt_.segment(dimv_, dim_kkt-dimv_) 
      = kkt_matrix_inverse_.block(dimx, dim_kkt-dimx, dim_kkt-dimx, dimx) 
          * x_res_;
  mu_tmp_active.noalias() -= dkkt_.segment(kkt_composition_.C_begin(), 
                                           kkt_composition_.C_size());
  a_tmp_.noalias() -= dkkt_.segment(kkt_composition_.Qa_begin(), 
                                    kkt_composition_.Qa_size());
  f_tmp_active.noalias() -= dkkt_.segment(kkt_composition_.Qf_begin(), 
                                          kkt_composition_.Qf_size());
  robot.integrateConfiguration(
      dkkt_.segment(kkt_composition_.Qq_begin(), kkt_composition_.Qq_size()), 
      -1, q_tmp_);
  v_tmp_.noalias() -= dkkt_.segment(kkt_composition_.Qv_begin(), 
                                    kkt_composition_.Qv_size());
}


void SplitParNMPC::forwardCollectionSerial(const Robot& robot,
                                           const Eigen::VectorXd& q_prev_old,   
                                           const Eigen::VectorXd& v_prev_old, 
                                           const Eigen::VectorXd& q_prev, 
                                           const Eigen::VectorXd& v_prev,
                                           Eigen::VectorXd& q, 
                                           Eigen::VectorXd& v) {
  robot.subtractConfiguration(q_prev, q_prev_old, x_res_.head(robot.dimv()));
  x_res_.tail(robot.dimv()) = v_prev - v_prev_old;
  dx_ = kkt_matrix_inverse_.block(kkt_composition_.Qx_begin(), 
                                  kkt_composition_.Fx_begin(), 
                                  kkt_composition_.Qx_size(), 
                                  kkt_composition_.Fx_size()) * x_res_;
  robot.integrateConfiguration(dx_.head(robot.dimv()), -1, q);
  v.noalias() -= dx_.tail(robot.dimv());
}


void SplitParNMPC::forwardCollectionParallel(const Robot& robot) {
  const int dim_kkt = kkt_matrix_.dimKKT();
  const int dimx = 2*robot.dimv();
  dkkt_.head(dim_kkt) 
      = kkt_matrix_inverse_.block(0, dim_kkt-dimx, dim_kkt, dimx) * x_res_;
  lmd_tmp_.noalias() -= dkkt_.segment(kkt_composition_.Fq_begin(), 
                                      kkt_composition_.Fq_size());
  gmm_tmp_.noalias() -= dkkt_.segment(kkt_composition_.Fv_begin(), 
                                      kkt_composition_.Fv_size());
  mu_tmp_active.noalias() -= dkkt_.segment(kkt_composition_.C_begin(), 
                                           kkt_composition_.C_size());
  a_tmp_.noalias() -= dkkt_.segment(kkt_composition_.Qa_begin(), 
                                    kkt_composition_.Qa_size());
  f_tmp_active.noalias() -= dkkt_.segment(kkt_composition_.Qf_begin(), 
                                          kkt_composition_.Qf_size());
  robot.integrateConfiguration(
      dkkt_.segment(kkt_composition_.Qq_begin(), kkt_composition_.Qq_size()), 
      -1, q_tmp_);
  v_tmp_.noalias() -= dkkt_.segment(kkt_composition_.Qv_begin(), 
                                    kkt_composition_.Qv_size());
}


void SplitParNMPC::computeDirection(const Robot& robot, 
                                    const Eigen::VectorXd& lmd, 
                                    const Eigen::VectorXd& gmm, 
                                    const Eigen::VectorXd& mu, 
                                    const Eigen::VectorXd& a,
                                    const Eigen::VectorXd& f, 
                                    const Eigen::VectorXd& q, 
                                    const Eigen::VectorXd& v, 
                                    const Eigen::VectorXd& u) {
  dlmd_ = lmd_tmp_ - lmd;
  dgmm_ = gmm_tmp_ - gmm;
  const int dimc = robot.dim_passive() + robot.dimf();
  dmu_.head(dimc) = mu_tmp_.head(dimc) - mu.head(dimc);
  da_ = a_tmp_ - a;
  df_.head(robot.dimf()) = f_tmp_.head(robot.dimf()) - f.head(robot.dimf());
  robot.subtractConfiguration(q_tmp_, q, dq_);
  dv_ = v_tmp_ - v;
  du_ = u_res_;
  du_.noalias() += du_dq_ * dq;
  du_.noalias() += du_dv_ * dv;
  du_.noalias() += du_da_ * da_;
  if (robot.dimf() > 0) {
    du_.noalias() += du_df_.leftCols(robot.dimf()) * df_.head(robot.dimf());
  }
  joint_constraints_.computeSlackAndDualDirection(dtau, dq, dv, da_, du_);
}

 
double SplitParNMPC::maxPrimalStepSize() {
  return joint_constraints_.maxSlackStepSize();
}


double SplitParNMPC::maxDualStepSize() {
  return joint_constraints_.maxDualStepSize();
}


std::pair<double, double> SplitParNMPC::costAndConstraintsViolation(
    Robot& robot, const double t, const double dtau, const Eigen::VectorXd& q, 
    const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
    const Eigen::VectorXd& f, const Eigen::VectorXd& u) {
  assert(dtau > 0);
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(u.size() == dimv_);
  assert(f.size() == max_dimf_);
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, q, v, a, f, u);
  cost += joint_constraints_.costSlackBarrier();
  double constraints_violation = 0;
  constraints_violation += q_res_.lpNorm<1>();
  constraints_violation += v_res_.lpNorm<1>();
  constraints_violation += dtau * u_res_.lpNorm<1>();
  constraints_violation += joint_constraints_.residualL1Nrom(dtau, q, v, a, u);
  constraints_violation += C_res_.head(dimc_).lpNorm<1>();
  return std::make_pair(cost, constraints_violation);
}


std::pair<double, double> SplitParNMPC::costAndConstraintsViolation(
    Robot& robot, const double step_size, const double t, const double dtau, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
    const Eigen::VectorXd& a, const Eigen::VectorXd& f, 
    const Eigen::VectorXd& u, const Eigen::VectorXd& q_next, 
    const Eigen::VectorXd& v_next, const Eigen::VectorXd& dq, 
    const Eigen::VectorXd& dv, const Eigen::VectorXd& dq_next, 
    const Eigen::VectorXd& dv_next) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(f.size() == max_dimf_);
  assert(u.size() == dimv_);
  assert(q_next.size() == dimq_);
  assert(v_next.size() == dimv_);
  assert(dq.size() == dimv_);
  assert(dv.size() == dimv_);
  assert(dq_next.size() == dimv_);
  assert(dv_next.size() == dimv_);
  q_tmp_ = q;
  robot.integrateConfiguration(dq, step_size, q_tmp_);
  v_tmp_ = v + step_size * dv;
  a_tmp_ = a + step_size * da_;
  u_tmp_ = u + step_size * du_;
  if (dimf_ > 0) {
    f_tmp_.head(dimf_) = f.head(dimf_) + step_size * df_.head(dimf_);
    robot.setContactForces(f_tmp_);
  }
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, q_tmp_, v_tmp_, a_tmp_, u_tmp_, f_tmp_);
  cost += joint_constraints_.costSlackBarrier(step_size);
  robot.subtractConfiguration(q_tmp_, q_next, q_res_);
  q_res_.noalias() += dtau * v_tmp_ - step_size * dq_next;
  v_res_ = v_tmp_ + dtau * a_tmp_ - v_next - step_size * dv_next;
  robot.RNEA(q_tmp_, v_tmp_, a_tmp_, u_res_tmp_);
  u_res_tmp_.noalias() -= u_tmp_;
  C_res_.head(dim_passive_) = dtau * u_tmp_.head(dim_passive_);
  if (dimf_ > 0) {
    robot.updateKinematics(q_tmp_, v_tmp_, a_tmp_);
    robot.computeBaumgarteResidual(dim_passive_, dtau, C_res_);
  }
  double constraints_violation = 0;
  constraints_violation += q_res_.lpNorm<1>();
  constraints_violation += v_res_.lpNorm<1>();
  constraints_violation += dtau * u_res_tmp_.lpNorm<1>();
  constraints_violation += joint_constraints_.residualL1Nrom(dtau, q_tmp_, 
                                                             v_tmp_, a_tmp_, 
                                                             u_tmp_);
  constraints_violation += C_res_.head(dimc_).lpNorm<1>();
  return std::make_pair(cost, constraints_violation);
}


void SplitParNMPC::updateDual(const double step_size) {
  assert(step_size > 0);
  assert(step_size <= 1);
  joint_constraints_.updateDual(step_size);
}


void SplitParNMPC::updatePrimal(Robot& robot, const double step_size, 
                                const double dtau, const Eigen::MatrixXd& Pqq, 
                                const Eigen::MatrixXd& Pqv, 
                                const Eigen::MatrixXd& Pvq, 
                                const Eigen::MatrixXd& Pvv, 
                                const Eigen::VectorXd& sq, 
                                const Eigen::VectorXd& sv, 
                                const Eigen::VectorXd& dq, 
                                const Eigen::VectorXd& dv, Eigen::VectorXd& lmd, 
                                Eigen::VectorXd& gmm, Eigen::VectorXd& q, 
                                Eigen::VectorXd& v, Eigen::VectorXd& a, 
                                Eigen::VectorXd& u, Eigen::VectorXd& beta, 
                                Eigen::VectorXd& f, Eigen::VectorXd& mu) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  assert(Pqq.rows() == dimv_);
  assert(Pqq.cols() == dimv_);
  assert(Pqv.rows() == dimv_);
  assert(Pqv.cols() == dimv_);
  assert(Pvq.rows() == dimv_);
  assert(Pvq.cols() == dimv_);
  assert(Pvv.rows() == dimv_);
  assert(Pvv.cols() == dimv_);
  assert(sq.size() == dimv_);
  assert(sv.size() == dimv_);
  assert(dq.size() == dimv_);
  assert(dv.size() == dimv_);
  assert(lmd.size() == dimv_);
  assert(gmm.size() == dimv_);
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(u.size() == dimv_);
  assert(beta.size() == dimv_);
  assert(f.size() == max_dimf_);
  assert(mu.size() == max_dimc_);
  lmd.noalias() += step_size * (Pqq * dq + Pqv * dv - sq);
  gmm.noalias() += step_size * (Pvq * dq + Pvv * dv - sv);
  robot.integrateConfiguration(dq, step_size, q);
  v.noalias() += step_size * dv;
  a.noalias() += step_size * da_;
  u.noalias() += step_size * du_;
  lu_.noalias() -= dtau * beta;
  beta.noalias() += step_size * lu_ / dtau;
  beta.noalias() += step_size * luu_ * du_ / dtau;
  if (dimf_ > 0) {
    f.head(dimf_).noalias() += step_size * df_.head(dimf_);
  }
  if (dimc_ > 0) {
    mu.head(dimc_).noalias() += step_size * dmu_.head(dimc_);
  }
  joint_constraints_.updateSlack(step_size);
}


void SplitParNMPC::getStateFeedbackGain(Eigen::MatrixXd& Kq, 
                                        Eigen::MatrixXd& Kv) const {
  assert(Kq.cols() == dimv_);
  assert(Kq.rows() == dimv_);
  assert(Kv.cols() == dimv_);
  assert(Kv.rows() == dimv_);
  Kq = du_dq_ + du_da_ * Kaq_ + du_df_.leftCols(dimf_) * Kfq_.topRows(dimf_);
  Kv = du_dv_ + du_da_ * Kav_ + du_df_.leftCols(dimf_) * Kfv_.topRows(dimf_);
}


double SplitParNMPC::squaredKKTErrorNorm(Robot& robot, const double t, 
                                         const double dtau, 
                                         const Eigen::VectorXd& q_prev, 
                                         const Eigen::VectorXd& v_prev, 
                                         const Eigen::VectorXd& lmd, 
                                         const Eigen::VectorXd& gmm, 
                                         const Eigen::VectorXd& mu, 
                                         const Eigen::VectorXd& a, 
                                         const Eigen::VectorXd& f, 
                                         const Eigen::VectorXd& q, 
                                         const Eigen::VectorXd& v, 
                                         const Eigen::VectorXd& u, 
                                         const Eigen::VectorXd& lmd_next,
                                         const Eigen::VectorXd& gmm_next,
                                         const Eigen::VectorXd& q_next,
                                         const Eigen::VectorXd& v_next) {
  assert(dtau > 0);
  assert(q_prev.size() == dimq_);
  assert(v_prev.size() == dimv_);
  assert(lmd.size() == dimv_);
  assert(gmm.size() == dimv_);
  assert(mu.size() == max_dimc_);
  assert(a.size() == dimv_);
  assert(f.size() == max_dimf_);
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(u.size() == dimv_);
  assert(lmd_next.size() == dimv_);
  assert(gmm_next.size() == dimv_);
  assert(q_next.size() == dimq_);
  assert(v_next.size() == dimv_);
  // Reset the KKT matrix and KKT residual.
  kkt_matrix_.setZero();
  kkt_matrix_.setContactStatus(robot);
  kkt_residual_.setZero();
  kkt_residual_.setContactStatus(robot);
  // Extract active blocks for matrices with associated with contacts.
  const Eigen::Ref<const Eigen::MatrixXd> du_df_active 
      = du_df_.leftCols(robot.dimf());
  const Eigen::Ref<const Eigen::VectorXd> mu_active 
      = mu.head(robot.dim_passive()+robot.dimf());
  const Eigen::Ref<const Eigen::VectorXd> mu_passive 
      = mu.head(robot.dim_passive());
  const Eigen::Ref<const Eigen::VectorXd> f_active = f.head(robot.dimf());
  Eigen::Ref<Eigen::VectorXd> mu_tmp_acitve = mu_tmp_.head(robot.dim_passive());
  Eigen::Ref<Eigen::VectorXd> f_tmp_active = f_tmp_.head(robot.dimf());
  // Compute KKT residual and KKT matrix.
  if (robot.dimf() > 0) {
    robot.updateKinematics(q, v, a);
  }
  parnmpclinearizer::linearizeStageCost(robot, cost_, cost_data_, t, dtau, 
                                        q, v, a, f, u, kkt_residual_, lu_);
  parnmpclinearizer::linearizeDynamics(robot, dtau, q, v, a, f, u,  
                                       q_prev, v_prev, kkt_residual_, u_res_, 
                                       du_dq_, du_dv_, du_da_, du_df_);
  parnmpclinearizer::linearizeConstraints(robot, dtau, u, u_res_,  
                                          du_dq_, du_dv_, du_da_, du_df_, 
                                          kkt_residual_, kkt_matrix_);
  // Augmnet the partial derivatives of the state equation.
  if (robot.has_floating_base()) {
    robot.dSubtractdConfigurationMinus(q_prev, q, dsubtract_dqminus_);
    robot.dSubtractdConfigurationPlus(q, q_next, dsubtract_dqplus_);
    kkt_residual_.lq().noalias() 
        += dsubtract_dqplus_.transpose() * lmd_next
            + dsubtract_dqminus_.transpose() * lmd;
    kkt_residual_.lv().noalias() += dtau * lmd - gmm + gmm_next;
    kkt_residual_.la().noalias() += dtau * gmm;
  }
  else {
    kkt_residual_.lq().noalias() += lmd_next - lmd;
    kkt_residual_.lv().noalias() += dtau * lmd - gmm + gmm_next;
    kkt_residual_.la().noalias() += dtau * gmm;
  }
  // Augmnet the partial derivatives of the inequality constriants.
  joint_constraints_.augmentDualResidual(dtau, kkt_residual_.lq(), 
                                         kkt_residual_.lv(), 
                                         kkt_residual_.la());
  joint_constraints_.augmentDualResidual(dtau, lu_);
  // Augment the equality constraints 
  kkt_residual_.lq().noalias() += kkt_matrix_.Cq().transpose() * mu_active;
  kkt_residual_.lv().noalias() += kkt_matrix_.Cv().transpose() * mu_active;
  kkt_residual_.la().noalias() += kkt_matrix_.Ca().transpose() * mu_active;
  if (dimf_ > 0) {
    kkt_residual_.lf().noalias() += kkt_matrix_.Cf().transpose() * mu_passive;
  }
  double error = 0;
  error += kkt_residual_.KKT_residual().squaredNorm();
  error += lu_.squaredNorm();
  error += joint_constraints_.residualSquaredNrom(dtau, q, v, a, u);
  return error;
}

} // namespace idocp