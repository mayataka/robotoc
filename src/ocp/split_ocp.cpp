#include "idocp/ocp/split_ocp.hpp"

#include <assert.h>


namespace idocp {

SplitOCP::SplitOCP(const Robot& robot, 
                   const std::shared_ptr<CostFunction>& cost,
                   const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(robot),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot)),
    kkt_residual_(robot),
    kkt_matrix_(robot),
    riccati_matrix_factorizer_(robot),
    riccati_matrix_inverter_(robot),
    s_tmp_(robot),
    ka_(Eigen::VectorXd::Zero(robot.dimv())),
    kf_(Eigen::VectorXd::Zero(robot.max_dimf())),
    kmu_(Eigen::VectorXd::Zero(robot.max_dimf()+robot.dim_passive())),
    Kaq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Kav_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Kfq_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
    Kfv_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
    Kmuq_(Eigen::MatrixXd::Zero(robot.max_dimf()+robot.dim_passive(),
                                robot.dimv())),
    Kmuv_(Eigen::MatrixXd::Zero(robot.max_dimf()+robot.dim_passive(), 
                                robot.dimv())) {
}


SplitOCP::SplitOCP() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    kkt_residual_(),
    kkt_matrix_(),
    riccati_matrix_factorizer_(),
    riccati_matrix_inverter_(),
    s_tmp_(),
    ka_(),
    kf_(),
    kmu_(),
    Kaq_(),
    Kav_(),
    Kfq_(),
    Kfv_(),
    Kmuq_(),
    Kmuv_() {
}


SplitOCP::~SplitOCP() {
}


bool SplitOCP::isFeasible(const Robot& robot, const SplitSolution& s) {
  return constraints_->isFeasible(robot, constraints_data_, s);
}


void SplitOCP::initConstraints(const Robot& robot, const int time_step, 
                               const double dtau, const SplitSolution& s) { 
  assert(time_step >= 0);
  assert(dtau > 0);
  constraints_->setSlackAndDual(robot, constraints_data_, dtau, s);
}


void SplitOCP::linearizeOCP(Robot& robot, const double t, const double dtau, 
                            const SplitSolution& s, 
                            const SplitSolution& s_next) {
  assert(dtau > 0);
  dimf_ = robot.dimf();
  dimc_ = robot.dim_passive() + robot.dimf();
  if (dimf_ > 0) {
    robot.updateKinematics(q, v, a);
  }
  kkt_residual_.setContactStatus(robot);
  kkt_residual_.setZeroMinimum();
  kkt_matrix_.setContactStatus(robot);
  kkt_matrix_.setZeroMinimum();
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, kkt_residual_);
  state_equation_.linearizeForwardEuler(robot, dtau, s, s_next, kkt_residual_);
  equalityconstraints::LinearizeEqualityConstraints(robot, dtau, s, 
                                                    kkt_matrix_, kkt_residual_);
  inverse_dynamics_.linearizeInverseDynamics(robot, dtau, s, kkt_residual_);

  joint_constraints_.augmentDualResidual(dtau, lu_);
  cost_->luu(robot, cost_data_, t, dtau, u, luu_);
  joint_constraints_.condenseSlackAndDual(dtau, u, luu_, lu_);
  lu_condensed_ = lu_ + luu_ * u_res_;
  // Augment the condensed Newton residual of the contorl input torques. 
  lq_.noalias() += du_dq_.transpose() * lu_condensed_;
  lv_.noalias() += du_dv_.transpose() * lu_condensed_;
  la_.noalias() += du_da_.transpose() * lu_condensed_;
  if (dimf_ > 0) {
    lf_.head(dimf_).noalias() 
                  += du_df_.leftCols(dimf_).transpose() * lu_condensed_;
  }
  // Augmnet the partial derivatives of the state equation.
  lq_.noalias() += lmd_next - lmd;
  lv_.noalias() += dtau * lmd_next + gmm_next - gmm;
  la_.noalias() += dtau * gmm_next;
  // Augmnet the partial derivatives of the inequality constriants.
  joint_constraints_.augmentDualResidual(dtau, lq_, lv_, la_);
  // Augment the equality constraints 
  lq_.noalias() += Cq_.topRows(dimc_).transpose() * mu.head(dimc_);
  lv_.noalias() += Cv_.topRows(dimc_).transpose() * mu.head(dimc_);
  la_.noalias() += Ca_.topRows(dimc_).transpose() * mu.head(dimc_);
  if (dimf_ > 0) {
    lf_.head(dimf_).noalias() += Cf_.leftCols(dimf_).transpose() 
        * mu.head(dim_passive_);
  }
  // Augment the condensed Hessian of the contorl input torques. 
  Qqq_ = du_dq_.transpose() * luu_ * du_dq_;
  Qqv_ = du_dq_.transpose() * luu_ * du_dv_;
  Qqa_ = du_dq_.transpose() * luu_ * du_da_;
  Qvv_ = du_dv_.transpose() * luu_ * du_dv_;
  Qva_ = du_dv_.transpose() * luu_ * du_da_;
  Qaa_ = du_da_.transpose() * luu_ * du_da_;
  Qvq_ = Qqv_.transpose();
  if (dimf_ > 0) {
    Qqf_.leftCols(dimf_) = du_dq_.transpose() * luu_ * du_df_.leftCols(dimf_);
    Qvf_.leftCols(dimf_) = du_dv_.transpose() * luu_ * du_df_.leftCols(dimf_);
    Qaf_.leftCols(dimf_) = du_da_.transpose() * luu_ * du_df_.leftCols(dimf_);
    Qff_.topLeftCorner(dimf_, dimf_) 
        = du_df_.leftCols(dimf_).transpose() * luu_ * du_df_.leftCols(dimf_);
  }
  // Condense the slack and dual variables of the inequality constraints on 
  // the configuration, velocity, and acceleration.
  joint_constraints_.condenseSlackAndDual(dtau, q, v, a, Qqq_, Qvv_, Qaa_, 
                                          lq_, lv_, la_);
  // Augment the cost function Hessian. 
  cost_->augment_lqq(robot, cost_data_, t, dtau, q, v, a, Qqq_);
  cost_->augment_lvv(robot, cost_data_, t, dtau, q, v, a, Qvv_);
  cost_->augment_laa(robot, cost_data_, t, dtau, q, v, a, Qaa_);
  if (dimf_ > 0) {
    cost_->augment_lff(robot, cost_data_, t, dtau, f, Qff_);
  }
  if (robot.has_floating_base()) {
    riccati_matrix_factorizer_.setIntegrationSensitivities(robot, dtau, q, v);
  }
  if (dimf_ > 0) {
    riccati_matrix_inverter_.setContactStatus(robot);
    riccati_matrix_inverter_.precompute(Qaf_, Qff_);
  }
}


void SplitOCP::backwardRiccatiRecursion(
    const double dtau, const RiccatiFactorization& riccati_next, 
    RiccatiFactorization& riccati) {
  assert(dtau > 0);
  // Qqq_, Qqv_, Qvq_, Qvv_: representing Riccati factorization F
  // Qqa_, Qqf_, Qva_, Qvf_ : representing Riccati factorization H
  // Qaa_, Qaf_, Qff_ : representing Riccati factorization G
  riccati_factorizer_.factorizeF(dtau, riccati_next.Pqq, riccati_next.Pqv, 
                                 riccati_next.Pvq, riccati_next.Pvv, 
                                 kkt_matrix_.Qqq(), kkt_matrix_.Qqv(), 
                                 kkt_matrix_.Qvq(), kkt_matrix_.Qvv());
  riccati_actorizer_.factorizeH(dtau, riccati_next.Pqv, riccati_next.Pvv, 
                                kkt_matrix_.Qqa(), kkt_matrix_.Qva());
  riccati_actorizer_.factorizeG(dtau, riccati_next.Pvv, kkt_matrix_.Qaa());
  kkt_residual_.la().noalias() += dtau * riccati_next.Pvq * kkt_residual_.Fq(); 
  kkt_residual_.la().noalias() += dtau * riccati_next.Pvv * kkt_residual_.Fv();
  kkt_residual_.la().noalias() -= dtau * riccati_next.sv;
  // Computes the state feedback gain and feedforward terms
  riccati_inverter_.invert(robot, kkt_matrix_.Qafaf(), kkt_matrix_.Caf());
  // Computes the Riccati factorization matrices
  riccati.Pqq = kkt_matrix_.Qqq();
  riccati.Pqq.noalias() += Kaq_.transpose() * kkt_matrix_.Qqa().transpose();
  riccati.Pqv = kkt_matrix_.Qqv();
  riccati.Pqv.noalias() += Kaq_.transpose() * kkt_matrix_.Qva().transpose();
  riccati.Pvq = kkt_matrix_.Qvq();
  riccati.Pvq.noalias() += Kav_.transpose() * kkt_matrix_.Qqa().transpose();
  riccati.Pvv = kkt_matrix_.Qvv();
  riccati.Pvv.noalias() += Kav_.transpose() * kkt_matrix_.Qva().transpose();
  // Computes the Riccati factorization vectors
  riccati.sq = riccati_next.sq - kkt_residual_.lq();
  riccati.sq.noalias() -= riccati_next.Pqq * kkt_residual_.Fq();
  riccati.sq.noalias() -= riccati_next.Pqv * kkt_residual_.Fv();
  riccati.sq.noalias() -= kkt_matrix_.Qqa() * ka_;
  riccati.sv = dtau * riccati_next.sq + riccati_next.sv - kkt_residual_.lv();
  riccati.sv.noalias() -= dtau * riccati_next.Pqq * kkt_residual_.Fq();
  riccati.sv.noalias() -= riccati_next.Pvq * kkt_residual_.Fq();
  riccati.sv.noalias() -= dtau * riccati_next.Pqv * kkt_residual_.Fv();
  riccati.sv.noalias() -= riccati_next.Pvv * kkt_residual_.Fv();
  riccati.sv.noalias() -= kkt_matrix_.Qva() * ka_;
  if (dimf_ > 0) {
    riccati.Pqq.noalias() += Kfq_active().transpose() * kkt_matrix_.Qqf().transpose();
    riccati.Pqv.noalias() += Kfq_active().transpose() * kkt_matrix_.Qvf().transpose();
    riccati.Pvq.noalias() += Kfv_active().transpose() * kkt_matrix_.Qqf().transpose();
    riccati.Pvv.noalias() += Kfv_active().transpose() * kkt_matrix_.Qvf().transpose();
    riccati.sq.noalias() -= kkt_matrix_.Qqf() * kf_active();
    riccati.sv.noalias() -= kkt_matrix_.Qvf() * kf_active();
  }
  if (dimc_ > 0) {
    riccati.Pqq.noalias() += Kmuq_active().transpose() * kkt_matrix_.Cq();
    riccati.Pqv.noalias() += Kmuq_active().transpose() * kkt_matrix_.Cv();
    riccati.Pvq.noalias() += Kmuv_active().transpose() * kkt_matrix_.Cq();
    riccati.Pvv.noalias() += Kmuv_active().transpose() * kkt_matrix_.Cv();
    riccati.sq.noalias() -= kkt_matrix_.Cq().transpose() * kmu_active();
    riccati.sv.noalias() -= kkt_matrix_.Cv().transpose() * kmu_active();
  }
}


void SplitOCP::forwardRiccatiRecursion(const double dtau, 
                                       const Eigen::VectorXd& dq,   
                                       const Eigen::VectorXd& dv, 
                                       Eigen::VectorXd& dq_next,
                                       Eigen::VectorXd& dv_next) {
  assert(dtau > 0);
  assert(dq.size() == dimv_);
  assert(dv.size() == dimv_);
  assert(dq_next.size() == dimv_);
  assert(dv_next.size() == dimv_);
  da_ = ka_ + Kaq_ * dq + Kav_ * dv;
  dq_next = dq + dtau * dv + q_res_;
  dv_next = dv + dtau * da_ + v_res_;
}


void SplitOCP::computeCondensedDirection(const double dtau, 
                                         const Eigen::VectorXd& dq, 
                                         const Eigen::VectorXd& dv) {
  assert(dtau > 0);
  assert(dq.size() == dimv_);
  assert(dv.size() == dimv_);
  if (dimf_ > 0) {
    df_.head(dimf_) = kf_.head(dimf_) + Kfq_.topRows(dimf_) * dq 
                                      + Kfv_.topRows(dimf_) * dv;
  }
  if (dimc_ > 0) {
    dmu_.head(dimc_) = kmu_.head(dimc_) + Kmuq_.topRows(dimc_) * dq 
                                        + Kmuv_.topRows(dimc_) * dv;
  }
  du_ = u_res_;
  du_.noalias() += du_dq_ * dq;
  du_.noalias() += du_dv_ * dv;
  du_.noalias() += du_da_ * da_;
  if (dimf_ > 0) {
    du_.noalias() += du_df_.leftCols(dimf_) * df_.head(dimf_);
  }
  joint_constraints_.computeSlackAndDualDirection(dtau, dq, dv, da_, du_);
}

 
double SplitOCP::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


double SplitOCP::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


std::pair<double, double> SplitOCP::costAndConstraintsViolation(
    Robot& robot, const double t, const double dtau, const Eigen::VectorXd& q, 
    const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
    const Eigen::VectorXd& u, const Eigen::VectorXd& f) {
  assert(dtau > 0);
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(u.size() == dimv_);
  assert(f.size() == max_dimf_);
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, q, v, a, u, f);
  cost += joint_constraints_.costSlackBarrier();
  double constraints_violation = 0;
  constraints_violation += q_res_.lpNorm<1>();
  constraints_violation += v_res_.lpNorm<1>();
  constraints_violation += dtau * u_res_.lpNorm<1>();
  constraints_violation += joint_constraints_.residualL1Nrom(dtau, q, v, a, u);
  constraints_violation += C_res_.head(dimc_).lpNorm<1>();
  return std::make_pair(cost, constraints_violation);
}


std::pair<double, double> SplitOCP::costAndConstraintsViolation(
    Robot& robot, const double step_size, const double t, const double dtau, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
    const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
    const Eigen::VectorXd& f, const Eigen::VectorXd& q_next, 
    const Eigen::VectorXd& v_next, const Eigen::VectorXd& dq, 
    const Eigen::VectorXd& dv, const Eigen::VectorXd& dq_next, 
    const Eigen::VectorXd& dv_next) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(u.size() == dimv_);
  assert(f.size() == max_dimf_);
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


void SplitOCP::updateDual(const double step_size) {
  assert(step_size > 0);
  assert(step_size <= 1);
  constraints_->updateDual(constraints_data_, step_size);
}


void SplitOCP::updatePrimal(Robot& robot, const double step_size, 
                            const double dtau,  
                            const RiccatiFactorization& riccati,
                            const Eigen::VectorXd& dq, 
                            const Eigen::VectorXd& dv, SplitSolution& s) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  assert(dq.size() == dimv_);
  assert(dv.size() == dimv_);
  s.lmd.noalias() 
      += step_size * (riccati.Pqq * dq + riccati.Pqv * dv - riccati.sq);
  s.gmm.noalias() 
      += step_size * (riccati.Pvq * dq + riccati.Pvv * dv - riccati.sv);
  robot.integrateConfiguration(dq, step_size, s.q);
  s.v.noalias() += step_size * dv;
  s.a.noalias() += step_size * da_;
  s.u.noalias() += step_size * du_;
  kkt_residual_.lu.noalias() -= dtau * s.beta;
  beta.noalias() += step_size * kkt_residual_.lu / dtau;
  beta.noalias() += step_size * luu_ * du_ / dtau;
  if (dimf_ > 0) {
    f.head(dimf_).noalias() += step_size * df_.head(dimf_);
  }
  if (dimc_ > 0) {
    mu.head(dimc_).noalias() += step_size * dmu_.head(dimc_);
  }
  joint_constraints_.updateSlack(step_size);
}


void SplitOCP::getStateFeedbackGain(Eigen::MatrixXd& Kq, 
                                    Eigen::MatrixXd& Kv) const {
  assert(Kq.cols() == dimv_);
  assert(Kq.rows() == dimv_);
  assert(Kv.cols() == dimv_);
  assert(Kv.rows() == dimv_);
  Kq = du_dq_ + du_da_ * Kaq_ + du_df_.leftCols(dimf_) * Kfq_.topRows(dimf_);
  Kv = du_dv_ + du_da_ * Kav_ + du_df_.leftCols(dimf_) * Kfv_.topRows(dimf_);
}


double SplitOCP::squaredKKTErrorNorm(Robot& robot, const double t, 
                                     const double dtau, const SplitSolution& s,
                                     const SplitSolution& s_next) {
  assert(dtau > 0);
  robot.updateKinematics(s.q, s.v, s.a);
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, 
                                    kkt_residual_);
  state_equation_.linearizeForwardEuler(robot, dtau, s, s_next, kkt_residual_);
  equalityconstraints::LinearizeEqualityConstraints(robot, dtau, s, 
                                                    kkt_matrix_, kkt_residual_);
  inverse_dynamics_.linearizeInverseDynamics(robot, dtau, s, kkt_residual_);
  double error = kkt_residual_.squaredKKTErrorNorm(dtau);
  error += constraints_->squaredKKTErrorNorm(robot, constraints_data_, dtau, s);
  return error;
}

} // namespace idocp