#include "ocp/split_ocp.hpp"

#include "Eigen/LU"


namespace invdynocp {

SplitOCP::SplitOCP(const Robot& robot, const CostFunctionInterface* cost, 
                   const ConstraintsInterface* constraints) 
  : cost_(const_cast<CostFunctionInterface*>(cost)),
    constraints_(const_cast<ConstraintsInterface*>(constraints)),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    u_(Eigen::VectorXd::Zero(robot.dimv())),
    lu_(Eigen::VectorXd::Zero(robot.dimv())),
    lq_(Eigen::VectorXd::Zero(robot.dimv())),
    lv_(Eigen::VectorXd::Zero(robot.dimv())),
    la_(Eigen::VectorXd::Zero(robot.dimv())),
    q_res_(Eigen::VectorXd::Zero(robot.dimv())),
    v_res_(Eigen::VectorXd::Zero(robot.dimv())),
    a_res_(Eigen::VectorXd::Zero(robot.dimv())),
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


void SplitOCP::linearizeOCP(Robot& robot, const double t, const double dtau, 
                            const Eigen::VectorXd& lmd, 
                            const Eigen::VectorXd& gmm, 
                            const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                            const Eigen::VectorXd& a, 
                            const Eigen::VectorXd& lmd_next, 
                            const Eigen::VectorXd& gmm_next, 
                            const Eigen::VectorXd& q_next,
                            const Eigen::VectorXd& v_next) {
  robot.RNEA(q, v, a, u_);
  cost_->lu(&robot, t, dtau, q, v, a, u_, lu_);
  cost_->lq(&robot, t, dtau, q, v, a, u_, lq_);
  cost_->lv(&robot, t, dtau, q, v, a, u_, lv_);
  cost_->la(&robot, t, dtau, q, v, a, u_, la_);
  robot.RNEADerivatives(q, v, a, du_dq_, du_dv_, du_da_);
  lq_.noalias() += du_dq_.transpose() * lu_;
  lv_.noalias() += du_dv_.transpose() * lu_;
  la_.noalias() += du_da_.transpose() * lu_;
  lq_.noalias() += lmd_next - lmd;
  lv_.noalias() += dtau * lmd_next + gmm_next - gmm;
  la_.noalias() += dtau * gmm_next;
  q_res_ = q + dtau * v - q_next;
  v_res_ = v + dtau * a - v_next;
  cost_->luu(&robot, t, dtau, q, v, a, u_, luu_);
  Qqq_ = luu_ * du_dq_.transpose() * du_dq_;
  Qqv_ = luu_ * du_dq_.transpose() * du_dv_;
  Qqa_ = luu_ * du_dq_.transpose() * du_da_;
  Qvv_ = luu_ * du_dv_.transpose() * du_dv_;
  Qva_ = luu_ * du_dv_.transpose() * du_da_;
  Qaa_ = luu_ * du_da_.transpose() * du_da_;
  Qvq_ = Qqv_.transpose();
  cost_->lqq(&robot, t, dtau, q, v, a, u_, Qqq_);
  cost_->lvv(&robot, t, dtau, q, v, a, u_, Qvv_);
  cost_->laa(&robot, t, dtau, q, v, a, u_, Qaa_);
}


void SplitOCP::linearizeTerminalCost(Robot& robot, const double t, 
                                     const Eigen::VectorXd& q, 
                                     const Eigen::VectorXd& v, 
                                     Eigen::MatrixXd& Qqq, Eigen::MatrixXd& Qqv, 
                                     Eigen::MatrixXd& Qvq, Eigen::MatrixXd& Qvv, 
                                     Eigen::VectorXd& Qq, Eigen::VectorXd& Qv) {
  cost_->phiq(&robot, t, q, v, Qq);
  cost_->phiv(&robot, t, q, v, Qv);
  Qq *= -1;
  Qv *= -1;
  cost_->phiqq(&robot, t, q, v, Qqq);
  cost_->phivv(&robot, t, q, v, Qvv);
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
  Ginv_ = Qaa_.inverse();
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
  da = k_ + Kq_ * dq + Kv_ * dv;
  dq_next = dq + dtau * dv + q_res_;
  dv_next = dv + dtau * da + v_res_;
}


void SplitOCP::updateOCP(Robot& robot, const Eigen::VectorXd& dq, 
                         const Eigen::VectorXd& dv, const Eigen::VectorXd& da, 
                         const Eigen::MatrixXd& Pqq, const Eigen::MatrixXd& Pqv, 
                         const Eigen::MatrixXd& Pvq, const Eigen::MatrixXd& Pvv, 
                         const Eigen::VectorXd& sq, const Eigen::VectorXd& sv, 
                         Eigen::VectorXd& q, Eigen::VectorXd& v, 
                         Eigen::VectorXd& a, Eigen::VectorXd& lmd, 
                         Eigen::VectorXd& gmm) {
  q.noalias() += dq;
  v.noalias() += dv;
  a.noalias() += da;
  lmd.noalias() += Pqq * dq + Pqv * dv - sq;
  gmm.noalias() += Pvq * dq + Pvv * dv - sv;
}


void SplitOCP::updateOCP(Robot& robot, const Eigen::VectorXd& dq, 
                         const Eigen::VectorXd& dv, const Eigen::MatrixXd& Pqq, 
                         const Eigen::MatrixXd& Pqv, const Eigen::MatrixXd& Pvq, 
                         const Eigen::MatrixXd& Pvv, const Eigen::VectorXd& sq, 
                         const Eigen::VectorXd& sv, Eigen::VectorXd& q, 
                         Eigen::VectorXd& v, Eigen::VectorXd& lmd, 
                         Eigen::VectorXd& gmm) {
  q.noalias() += dq;
  v.noalias() += dv;
  lmd.noalias() += Pqq * dq + Pqv * dv - sq;
  gmm.noalias() += Pvq * dq + Pvv * dv - sv;
}


double SplitOCP::optimalityError(Robot& robot, const double t, 
                                 const double dtau, const Eigen::VectorXd& lmd, 
                                 const Eigen::VectorXd& gmm, 
                                 const Eigen::VectorXd& q, 
                                 const Eigen::VectorXd& v, 
                                 const Eigen::VectorXd& a, 
                                 const Eigen::VectorXd& lmd_next, 
                                 const Eigen::VectorXd& gmm_next, 
                                 const Eigen::VectorXd& q_next,
                                 const Eigen::VectorXd& v_next) {
  robot.RNEA(q, v, a, u_);
  cost_->lu(&robot, t, dtau, q, v, a, u_, lu_);
  cost_->lq(&robot, t, dtau, q, v, a, u_, lq_);
  cost_->lv(&robot, t, dtau, q, v, a, u_, lv_);
  cost_->la(&robot, t, dtau, q, v, a, u_, la_);
  robot.RNEADerivatives(q, v, a, du_dq_, du_dv_, du_da_);
  lq_.noalias() += du_dq_.transpose() * lu_;
  lv_.noalias() += du_dv_.transpose() * lu_;
  la_.noalias() += du_da_.transpose() * lu_;
  lq_.noalias() += lmd_next - lmd;
  lv_.noalias() += dtau * lmd_next + gmm_next - gmm;
  la_.noalias() += dtau * gmm_next;
  q_res_ = q + dtau * v - q_next;
  v_res_ = v + dtau * a - v_next;
  double error = 0;
  error += lq_.norm();
  error += lv_.norm();
  error += la_.norm();
  error += q_res_.norm();
  error += v_res_.norm();
  return error;
}


double SplitOCP::terminalError(Robot& robot, const double t, 
                               const Eigen::VectorXd& lmd, 
                               const Eigen::VectorXd& gmm, 
                               const Eigen::VectorXd& q, 
                               const Eigen::VectorXd& v) {
  cost_->phiq(&robot, t, q, v, lq_);
  cost_->phiv(&robot, t, q, v, lv_);
  lq_.noalias() -= lmd;
  lv_.noalias() -= gmm;
  double error = lq_.norm() + lv_.norm();
  return error;
}


} // namespace invdynocp