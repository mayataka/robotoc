#include "ocp/split_OCP.hpp"

#include <assert.h>
#include "Eigen/LU"


namespace idocp {

SplitOCP::SplitOCP(const Robot& robot, const CostFunctionInterface* cost, 
                   const ConstraintsInterface* constraints) 
  : cost_(const_cast<CostFunctionInterface*>(cost)),
    constraints_(const_cast<ConstraintsInterface*>(constraints)),
    // joint_constraints_barrier_(robot),
    joint_constraints_(robot),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    lu_(Eigen::VectorXd::Zero(robot.dimv())),
    lq_(Eigen::VectorXd::Zero(robot.dimv())),
    lv_(Eigen::VectorXd::Zero(robot.dimv())),
    la_(Eigen::VectorXd::Zero(robot.dimv())),
    q_res_(Eigen::VectorXd::Zero(robot.dimv())),
    v_res_(Eigen::VectorXd::Zero(robot.dimv())),
    a_res_(Eigen::VectorXd::Zero(robot.dimv())),
    u_res_(Eigen::VectorXd::Zero(robot.dimv())),
    du_(Eigen::VectorXd::Zero(robot.dimv())),
    du_dq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    du_dv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    du_da_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    luu_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qqq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qqv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qqa_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qvq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qvv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qva_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qaa_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Ginv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Kq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Kv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    k_(Eigen::VectorXd::Zero(robot.dimv())) {
}


SplitOCP::~SplitOCP() {
}


bool SplitOCP::isFeasible(Robot& robot, const Eigen::VectorXd& q, 
                          const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                          const Eigen::VectorXd& u) {
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(a.size() == robot.dimv());
  assert(u.size() == robot.dimv());
  return joint_constraints_.isFeasible(robot, q, v, a, u);
}


void SplitOCP::initConstraints(Robot& robot, const double dtau, 
                               const Eigen::VectorXd& q, 
                               const Eigen::VectorXd& v, 
                               const Eigen::VectorXd& a, 
                               const Eigen::VectorXd& u) {
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(a.size() == robot.dimv());
  assert(u.size() == robot.dimv());
  joint_constraints_.setSlackAndDual(robot, dtau, q, v, a, u);
}


void SplitOCP::linearizeOCP(Robot& robot, const double t, const double dtau, 
                            const Eigen::VectorXd& lmd, 
                            const Eigen::VectorXd& gmm, 
                            const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                            const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                            const Eigen::VectorXd& beta, 
                            const Eigen::VectorXd& lmd_next, 
                            const Eigen::VectorXd& gmm_next, 
                            const Eigen::VectorXd& q_next,
                            const Eigen::VectorXd& v_next) {
  assert(dtau > 0);
  assert(lmd.size() == robot.dimv());
  assert(gmm.size() == robot.dimv());
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(a.size() == robot.dimv());
  assert(u.size() == robot.dimv());
  assert(beta.size() == robot.dimv());
  assert(lmd_next.size() == robot.dimv());
  assert(gmm_next.size() == robot.dimv());
  assert(q_next.size() == robot.dimq());
  assert(v_next.size() == robot.dimv());
  robot.RNEADerivatives(q, v, a, du_dq_, du_dv_, du_da_);
  cost_->lq(robot, t, dtau, q, v, a, lq_);
  lq_.noalias() += lmd_next - lmd;
  lq_.noalias() += dtau * du_dq_ * beta;
  cost_->lv(robot, t, dtau, q, v, a, lv_);
  lv_.noalias() += dtau * lmd_next + gmm_next - gmm;
  lv_.noalias() += dtau * du_dv_ * beta;
  cost_->la(robot, t, dtau, q, v, a, la_);
  la_.noalias() += dtau * gmm_next;
  la_.noalias() += dtau * du_da_ * beta;
  cost_->lu(robot, t, dtau, u, lu_);
  cost_->luu(robot, t, dtau, u, luu_);
  joint_constraints_.condenseSlackAndDual(robot, dtau, u, luu_, lu_);
  lu_.noalias() -= dtau * beta;
  robot.RNEA(q, v, a, u_res_);
  u_res_.noalias() -= u;
  du_ = lu_;
  du_.noalias() += luu_ * u_res_;
  lq_.noalias() += du_dq_.transpose() * du_;
  lv_.noalias() += du_dv_.transpose() * du_;
  la_.noalias() += du_da_.transpose() * du_;
  q_res_ = q + dtau * v - q_next;
  v_res_ = v + dtau * a - v_next;
  Qqq_ = du_dq_.transpose() * luu_ * du_dq_;
  Qqv_ = du_dq_.transpose() * luu_ * du_dv_;
  Qqa_ = du_dq_.transpose() * luu_ * du_da_;
  Qvv_ = du_dv_.transpose() * luu_ * du_dv_;
  Qva_ = du_dv_.transpose() * luu_ * du_da_;
  Qaa_ = du_da_.transpose() * luu_ * du_da_;
  joint_constraints_.condenseSlackAndDual(robot, dtau, q, v, a, Qqq_, Qvv_, 
                                          Qaa_, lq_, lv_, la_);
  cost_->lqq(robot, t, dtau, q, v, a, Qqq_);
  cost_->lvv(robot, t, dtau, q, v, a, Qvv_);
  cost_->laa(robot, t, dtau, q, v, a, Qaa_);
  Qvq_ = Qqv_.transpose();
}


void SplitOCP::linearizeOCP(Robot& robot, const double t, 
                            const Eigen::VectorXd& lmd, 
                            const Eigen::VectorXd& gmm, 
                            const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                            Eigen::MatrixXd& Qqq, Eigen::MatrixXd& Qqv, 
                            Eigen::MatrixXd& Qvq, Eigen::MatrixXd& Qvv, 
                            Eigen::VectorXd& Qq, Eigen::VectorXd& Qv) {
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(Qqq.rows() == robot.dimv());
  assert(Qqq.cols() == robot.dimv());
  assert(Qqv.rows() == robot.dimv());
  assert(Qqv.cols() == robot.dimv());
  assert(Qvq.rows() == robot.dimv());
  assert(Qvq.cols() == robot.dimv());
  assert(Qvv.rows() == robot.dimv());
  assert(Qvv.cols() == robot.dimv());
  assert(Qq.rows() == robot.dimv());
  assert(Qv.cols() == robot.dimv());
  cost_->phiq(robot, t, q, v, lq_);
  cost_->phiv(robot, t, q, v, lv_);
  Qq = lmd - lq_;
  Qv = gmm - lv_;
  cost_->phiqq(robot, t, q, v, Qqq);
  cost_->phivv(robot, t, q, v, Qvv);
}



void SplitOCP::backwardRecursion(const double dtau, 
                                 const Eigen::MatrixXd& Pqq_next, 
                                 const Eigen::MatrixXd& Pqv_next, 
                                 const Eigen::MatrixXd& Pvq_next, 
                                 const Eigen::MatrixXd& Pvv_next, 
                                 const Eigen::VectorXd& sq_next, 
                                 const Eigen::VectorXd& sv_next, 
                                 Eigen::MatrixXd& Pqq, Eigen::MatrixXd& Pqv, 
                                 Eigen::MatrixXd& Pvq, Eigen::MatrixXd& Pvv, 
                                 Eigen::VectorXd& sq, Eigen::VectorXd& sv) {
  assert(dtau > 0);
  assert(Pqq_next.rows() == Pqq_next.cols());
  assert(Pqv_next.rows() == Pqv_next.cols());
  assert(Pvq_next.rows() == Pvq_next.cols());
  assert(Pvv_next.rows() == Pvv_next.cols());
  assert(Pqq.rows() == Pqq.cols());
  assert(Pqv.rows() == Pqv.cols());
  assert(Pvq.rows() == Pvq.cols());
  assert(Pvv.rows() == Pvv.cols());
  assert(Pqq_next.rows() == Pqv_next.rows());
  assert(Pqv_next.rows() == Pvq_next.rows());
  assert(Pvq_next.rows() == Pvv_next.rows());
  assert(Pvv_next.rows() == sq_next.size());
  assert(sq_next.size() == sv_next.size());
  assert(sv_next.size() == Pqq.rows());
  assert(Pqq.rows() == Pqv.rows());
  assert(Pqv.rows() == Pvq.rows());
  assert(Pvq.rows() == Pvv.rows());
  assert(Pvv.rows() == sq.size());
  assert(sq.size() == sv.size());
  Qqq_.noalias() += Pqq_next;
  Qqv_.noalias() += dtau * Pqq_next;
  Qqv_.noalias() += Pqv_next;
  Qvq_.noalias() += dtau * Pqq_next;
  Qvq_.noalias() += Pvq_next;
  Qvv_.noalias() += (dtau*dtau) * Pqq_next;
  Qvv_.noalias() += dtau * (Pqv_next + Pvq_next);
  Qvv_.noalias() += Pvv_next;
  Qqa_.noalias() += dtau * Pqv_next;
  Qva_.noalias() += (dtau*dtau) * Pqv_next + dtau * Pvv_next;
  Qaa_.noalias() += (dtau*dtau) * Pvv_next;
  // Ginv_ = Qaa_.inverse();
  Ginv_ = Qaa_.llt().solve(Eigen::MatrixXd::Identity(dimv_, dimv_));
  Kq_ = - Ginv_ * Qqa_.transpose();
  Kv_ = - Ginv_ * Qva_.transpose();
  Pqq = Qqq_;
  Pqq.noalias() += Kq_.transpose() * Qqa_.transpose();
  Pqv = Qqv_;
  Pqv.noalias() += Kq_.transpose() * Qva_.transpose();
  Pvq = Qvq_;
  Pvq.noalias() += Kv_.transpose() * Qqa_.transpose();
  Pvv = Qvv_;
  Pvv.noalias() += Kv_.transpose() * Qva_.transpose();
  la_.noalias() += dtau * Pvq_next * q_res_;
  la_.noalias() += dtau * Pvv_next * v_res_;
  la_.noalias() -= dtau * sv_next;
  k_ = - Ginv_ * la_;
  sq = sq_next - lq_;
  sv = dtau * sq_next + sv_next - lv_;
  sq.noalias() -= Pqq_next * q_res_;
  sq.noalias() -= Pqv_next * v_res_;
  sv.noalias() -= dtau * Pqq_next * q_res_;
  sv.noalias() -= Pvq_next * q_res_;
  sv.noalias() -= dtau * Pqv_next * v_res_;
  sv.noalias() -= Pvv_next * v_res_;
  sq.noalias() -= Qqa_ * k_;
  sv.noalias() -= Qva_ * k_;
}


void SplitOCP::forwardRecursion(const double dtau, const Eigen::VectorXd& dq,   
                                const Eigen::VectorXd& dv, Eigen::VectorXd& da, 
                                Eigen::VectorXd& dq_next, 
                                Eigen::VectorXd& dv_next) const {
  assert(dtau > 0);
  assert(dq.size() == dv.size());
  assert(dv.size() == da.size());
  assert(da.size() == dq_next.size());
  assert(dq_next.size() == dv_next.size());
  da = k_ + Kq_ * dq + Kv_ * dv;
  dq_next = dq + dtau * dv + q_res_;
  dv_next = dv + dtau * da + v_res_;
}


std::pair<double, double> SplitOCP::computeMaxStepSize(
    Robot& robot, const double dtau, const Eigen::VectorXd& dq, 
    const Eigen::VectorXd& dv, const Eigen::VectorXd& da) {
  assert(dtau > 0);
  assert(dq.size() == robot.dimv());
  assert(dv.size() == robot.dimv());
  assert(da.size() == robot.dimv());
  du_ = u_res_;
  du_.noalias() += du_dq_ * dq;
  du_.noalias() += du_dv_ * dv;
  du_.noalias() += du_da_ * da;
  const std::pair<double, double> primal_dual_step_size 
      = joint_constraints_.computeDirectionAndMaxStepSize(robot, dtau, dq, dv, 
                                                          da, du_);
  return primal_dual_step_size;
}


std::pair<double, double> SplitOCP::computeCostAndConstraintsReisdual(
    Robot& robot, const double t, const double dtau, const Eigen::VectorXd& q, 
    const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
    const Eigen::VectorXd& u, const Eigen::VectorXd& q_next, 
    const Eigen::VectorXd& v_next) {
  assert(dtau > 0);
  assert(q.rows() == robot.dimq());
  assert(v.rows() == robot.dimv());
  assert(a.rows() == robot.dimv());
  assert(u.rows() == robot.dimv());
  assert(q_next.rows() == robot.dimq());
  assert(v_next.rows() == robot.dimv());
  double cost = 0;
  cost += cost_->l(robot, t, dtau, q, v, a, u);
  cost += joint_constraints_.slackBarrier();
  double constraints_residual = 0;
  q_res_ = q + dtau * v - q_next;
  v_res_ = v + dtau * a - v_next;
  constraints_residual += q_res_.lpNorm<1>();
  constraints_residual += v_res_.lpNorm<1>();
  constraints_residual += dtau * u_res_.lpNorm<1>();
  constraints_residual += joint_constraints_.residualL1Nrom(robot, dtau, q, v, 
                                                            a, u);
  return std::make_pair(cost, constraints_residual);
}


std::pair<double, double> SplitOCP::computeCostAndConstraintsReisdual(
    Robot& robot, const double step_size, const double t, const double dtau, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
    const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
    const Eigen::VectorXd& q_next, const Eigen::VectorXd& v_next, 
    const Eigen::VectorXd& dq, const Eigen::VectorXd& dv, 
    const Eigen::VectorXd& da, const Eigen::VectorXd& dq_next, 
    const Eigen::VectorXd& dv_next) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  assert(q.rows() == robot.dimq());
  assert(v.rows() == robot.dimv());
  assert(a.rows() == robot.dimv());
  assert(u.rows() == robot.dimv());
  assert(q_next.rows() == robot.dimq());
  assert(v_next.rows() == robot.dimv());
  assert(dq.size() == robot.dimv());
  assert(dv.size() == robot.dimv());
  assert(da.size() == robot.dimv());
  assert(dq_next.size() == robot.dimv());
  assert(dv_next.size() == robot.dimv());
  q_tmp_ = q + step_size * dq;
  v_tmp_ = v + step_size * dv;
  a_tmp_ = a + step_size * da;
  u_tmp_ = u + step_size * du_;
  double cost = 0;
  cost += cost_->l(robot, t, dtau, q_tmp_, v_tmp_, a_tmp_, u_tmp_);
  cost += joint_constraints_.slackBarrier(step_size);
  double constraints_residual = 0;
  q_res_ = q_tmp_ + dtau * v_tmp_ - q_next - step_size * dq_next;
  v_res_ = v_tmp_ + dtau * a_tmp_ - v_next - step_size * dv_next;
  constraints_residual += q_res_.lpNorm<1>();
  constraints_residual += v_res_.lpNorm<1>();
  constraints_residual += joint_constraints_.residualL1Nrom(robot, dtau, q_tmp_, 
                                                            v_tmp_, a_tmp_, 
                                                            u_tmp_);
  robot.RNEA(q_tmp_, v_tmp_, a_tmp_, u_res_);
  u_res_.noalias() -= u_tmp_;
  constraints_residual += dtau * u_res_.lpNorm<1>();
  return std::make_pair(cost, constraints_residual);
}


double SplitOCP::computeTerminalCost(Robot& robot, const double t, 
                                     const Eigen::VectorXd& q, 
                                     const Eigen::VectorXd& v) {
  return cost_->phi(robot, t, q, v);
}


double SplitOCP::computeTerminalCost(Robot& robot, const double step_size, 
                                     const double t, const Eigen::VectorXd& q, 
                                     const Eigen::VectorXd& v, 
                                     const Eigen::VectorXd& dq, 
                                     const Eigen::VectorXd& dv) {
  q_tmp_ = q + step_size * dq;
  v_tmp_ = v + step_size * dv;
  return cost_->phi(robot, t, q_tmp_, v_tmp_);
}


void SplitOCP::updateDual(const double step_size) {
  assert(step_size > 0);
  assert(step_size <= 1);
  joint_constraints_.updateDual(step_size);
}


void SplitOCP::updatePrimal(Robot& robot, const double step_size, 
                            const double dtau, const Eigen::VectorXd& dq, 
                            const Eigen::VectorXd& dv, 
                            const Eigen::VectorXd& da, 
                            const Eigen::MatrixXd& Pqq, 
                            const Eigen::MatrixXd& Pqv, 
                            const Eigen::MatrixXd& Pvq, 
                            const Eigen::MatrixXd& Pvv, 
                            const Eigen::VectorXd& sq, 
                            const Eigen::VectorXd& sv, Eigen::VectorXd& q, 
                            Eigen::VectorXd& v, Eigen::VectorXd& a, 
                            Eigen::VectorXd& u, Eigen::VectorXd& beta, 
                            Eigen::VectorXd& lmd, Eigen::VectorXd& gmm) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  assert(dq.size() == robot.dimv());
  assert(dv.size() == robot.dimv());
  assert(da.size() == robot.dimv());
  assert(Pqq.rows() == robot.dimv());
  assert(Pqq.cols() == robot.dimv());
  assert(Pqv.rows() == robot.dimv());
  assert(Pqv.cols() == robot.dimv());
  assert(Pvq.rows() == robot.dimv());
  assert(Pvq.cols() == robot.dimv());
  assert(Pvv.rows() == robot.dimv());
  assert(Pvv.cols() == robot.dimv());
  assert(sq.size() == robot.dimv());
  assert(sv.size() == robot.dimv());
  assert(q.rows() == robot.dimq());
  assert(v.rows() == robot.dimv());
  assert(a.rows() == robot.dimv());
  assert(lmd.rows() == robot.dimv());
  assert(gmm.rows() == robot.dimv());
  q.noalias() += step_size * dq;
  v.noalias() += step_size * dv;
  a.noalias() += step_size * da;
  u.noalias() += step_size * du_;
  beta.noalias() += step_size * lu_ / dtau;
  beta.noalias() += step_size * luu_ * du_ / dtau;
  lmd.noalias() += step_size * (Pqq * dq + Pqv * dv - sq);
  gmm.noalias() += step_size * (Pvq * dq + Pvv * dv - sv);
  joint_constraints_.updateSlack(step_size);
}


void SplitOCP::updatePrimal(Robot& robot, const double step_size, 
                            const Eigen::VectorXd& dq, 
                            const Eigen::VectorXd& dv, 
                            const Eigen::MatrixXd& Pqq, 
                            const Eigen::MatrixXd& Pqv, 
                            const Eigen::MatrixXd& Pvq, 
                            const Eigen::MatrixXd& Pvv, 
                            const Eigen::VectorXd& sq, 
                            const Eigen::VectorXd& sv, Eigen::VectorXd& q, 
                            Eigen::VectorXd& v, Eigen::VectorXd& lmd, 
                            Eigen::VectorXd& gmm) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dq.size() == robot.dimv());
  assert(dv.size() == robot.dimv());
  assert(Pqq.rows() == robot.dimv());
  assert(Pqq.cols() == robot.dimv());
  assert(Pqv.rows() == robot.dimv());
  assert(Pqv.cols() == robot.dimv());
  assert(Pvq.rows() == robot.dimv());
  assert(Pvq.cols() == robot.dimv());
  assert(Pvv.rows() == robot.dimv());
  assert(Pvv.cols() == robot.dimv());
  assert(sq.size() == robot.dimv());
  assert(sv.size() == robot.dimv());
  assert(q.rows() == robot.dimq());
  assert(v.rows() == robot.dimv());
  assert(lmd.rows() == robot.dimv());
  assert(gmm.rows() == robot.dimv());
  q.noalias() += step_size * dq;
  v.noalias() += step_size * dv;
  lmd.noalias() += step_size * (Pqq * dq + Pqv * dv - sq);
  gmm.noalias() += step_size * (Pvq * dq + Pvv * dv - sv);
}


void SplitOCP::updatePrimal(Robot& robot, const double dtau,   
                            const Eigen::VectorXd& dq, 
                            const Eigen::VectorXd& dv, 
                            const Eigen::VectorXd& da, 
                            const Eigen::MatrixXd& Pqq, 
                            const Eigen::MatrixXd& Pqv, 
                            const Eigen::MatrixXd& Pvq, 
                            const Eigen::MatrixXd& Pvv, 
                            const Eigen::VectorXd& sq, 
                            const Eigen::VectorXd& sv, Eigen::VectorXd& q, 
                            Eigen::VectorXd& v, Eigen::VectorXd& a, 
                            Eigen::VectorXd& u, Eigen::VectorXd& beta, 
                            Eigen::VectorXd& lmd, Eigen::VectorXd& gmm) {
  assert(dtau > 0);
  assert(dq.size() == robot.dimv());
  assert(dv.size() == robot.dimv());
  assert(da.size() == robot.dimv());
  assert(Pqq.rows() == robot.dimv());
  assert(Pqq.cols() == robot.dimv());
  assert(Pqv.rows() == robot.dimv());
  assert(Pqv.cols() == robot.dimv());
  assert(Pvq.rows() == robot.dimv());
  assert(Pvq.cols() == robot.dimv());
  assert(Pvv.rows() == robot.dimv());
  assert(Pvv.cols() == robot.dimv());
  assert(sq.size() == robot.dimv());
  assert(sv.size() == robot.dimv());
  assert(q.rows() == robot.dimq());
  assert(v.rows() == robot.dimv());
  assert(a.rows() == robot.dimv());
  assert(lmd.rows() == robot.dimv());
  assert(gmm.rows() == robot.dimv());
  q.noalias() += dq;
  v.noalias() += dv;
  a.noalias() += da;
  du_ = u_res_;
  du_.noalias() += du_dq_ * dq;
  du_.noalias() += du_dv_ * dv;
  du_.noalias() += du_da_ * da;
  u.noalias() += du_;
  beta.noalias() += lu_ / dtau;
  beta.noalias() += luu_ * du_ / dtau;
  lmd.noalias() += Pqq * dq + Pqv * dv - sq;
  gmm.noalias() += Pvq * dq + Pvv * dv - sv;
}


void SplitOCP::updatePrimal(Robot& robot, const Eigen::VectorXd& dq, 
                            const Eigen::VectorXd& dv, 
                            const Eigen::MatrixXd& Pqq, 
                            const Eigen::MatrixXd& Pqv, 
                            const Eigen::MatrixXd& Pvq, 
                            const Eigen::MatrixXd& Pvv, 
                            const Eigen::VectorXd& sq, 
                            const Eigen::VectorXd& sv, Eigen::VectorXd& q, 
                            Eigen::VectorXd& v, Eigen::VectorXd& lmd, 
                            Eigen::VectorXd& gmm) {
  assert(dq.size() == robot.dimv());
  assert(dv.size() == robot.dimv());
  assert(Pqq.rows() == robot.dimv());
  assert(Pqq.cols() == robot.dimv());
  assert(Pqv.rows() == robot.dimv());
  assert(Pqv.cols() == robot.dimv());
  assert(Pvq.rows() == robot.dimv());
  assert(Pvq.cols() == robot.dimv());
  assert(Pvv.rows() == robot.dimv());
  assert(Pvv.cols() == robot.dimv());
  assert(sq.size() == robot.dimv());
  assert(sv.size() == robot.dimv());
  assert(q.rows() == robot.dimq());
  assert(v.rows() == robot.dimv());
  assert(lmd.rows() == robot.dimv());
  assert(gmm.rows() == robot.dimv());
  q.noalias() += dq;
  v.noalias() += dv;
  lmd.noalias() += Pqq * dq + Pqv * dv - sq;
  gmm.noalias() += Pvq * dq + Pvv * dv - sv;
}


double SplitOCP::squaredOCPErrorNorm(Robot& robot, const double t, 
                                     const double dtau, 
                                     const Eigen::VectorXd& lmd, 
                                     const Eigen::VectorXd& gmm, 
                                     const Eigen::VectorXd& q, 
                                     const Eigen::VectorXd& v, 
                                     const Eigen::VectorXd& a, 
                                     const Eigen::VectorXd& u, 
                                     const Eigen::VectorXd& beta, 
                                     const Eigen::VectorXd& lmd_next, 
                                     const Eigen::VectorXd& gmm_next, 
                                     const Eigen::VectorXd& q_next,
                                     const Eigen::VectorXd& v_next) {
  assert(dtau > 0);
  assert(lmd.rows() == robot.dimv());
  assert(gmm.rows() == robot.dimv());
  assert(q.rows() == robot.dimq());
  assert(v.rows() == robot.dimv());
  assert(a.rows() == robot.dimv());
  assert(u.rows() == robot.dimv());
  assert(beta.rows() == robot.dimv());
  assert(lmd_next.rows() == robot.dimv());
  assert(gmm_next.rows() == robot.dimv());
  assert(q_next.rows() == robot.dimq());
  assert(v_next.rows() == robot.dimv());
  robot.RNEADerivatives(q, v, a, du_dq_, du_dv_, du_da_);
  cost_->lq(robot, t, dtau, q, v, a, lq_);
  lq_.noalias() += lmd_next - lmd;
  lq_.noalias() += dtau * du_dq_ * beta;
  // If there are equality constraints, augment to lq_

  cost_->lv(robot, t, dtau, q, v, a, lv_);
  lv_.noalias() += dtau * lmd_next + gmm_next - gmm;
  lv_.noalias() += dtau * du_dv_ * beta;
  // If there are equality constraints, augment to lv_

  cost_->la(robot, t, dtau, q, v, a, la_);
  la_.noalias() += dtau * gmm_next;
  la_.noalias() += dtau * du_da_ * beta;
  // If there are equality constraints, augment to la_

  cost_->lu(robot, t, dtau, u, lu_);
  lu_.noalias() -= dtau * beta;
  // If there are equality constraints, augment to lu_

  // std::cout << "dual error before = " << lq_.squaredNorm() + lv_.squaredNorm() + la_.squaredNorm() + lu_.squaredNorm() << std::endl;
  joint_constraints_.augmentDualResidual(robot, dtau, lq_, lv_, la_, lu_);
  // std::cout << "dual error after = " << lq_.squaredNorm() + lv_.squaredNorm() + la_.squaredNorm() + lu_.squaredNorm() << std::endl;

  q_res_ = q + dtau * v - q_next;
  v_res_ = v + dtau * a - v_next;
  robot.RNEA(q, v, a, u_res_);
  u_res_.noalias() -= u;
  u_res_.array() *= dtau;

  double error = 0;
  error += lq_.squaredNorm();
  error += lv_.squaredNorm();
  error += la_.squaredNorm();
  error += lu_.squaredNorm();
  error += q_res_.squaredNorm();
  error += v_res_.squaredNorm();
  error += u_res_.squaredNorm();
  error += joint_constraints_.residualSquaredNrom(robot, dtau, q, v, a, u);
  return error;
}


double SplitOCP::squaredOCPErrorNorm(Robot& robot, const double t, 
                                     const Eigen::VectorXd& lmd, 
                                     const Eigen::VectorXd& gmm, 
                                     const Eigen::VectorXd& q, 
                                     const Eigen::VectorXd& v) {
  assert(lmd.rows() == robot.dimv());
  assert(gmm.rows() == robot.dimv());
  assert(q.rows() == robot.dimq());
  assert(v.rows() == robot.dimv());
  cost_->phiq(robot, t, q, v, lq_);
  cost_->phiv(robot, t, q, v, lv_);
  lq_.noalias() -= lmd;
  lv_.noalias() -= gmm;
  double error = lq_.squaredNorm() + lv_.squaredNorm();
  return error;
}

} // namespace idocp