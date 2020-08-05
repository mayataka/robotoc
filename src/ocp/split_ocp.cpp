#include "idocp/ocp/split_ocp.hpp"
#include "idocp/ocp/ocp_linearizer.hpp"

#include <assert.h>


namespace idocp {

SplitOCP::SplitOCP(const Robot& robot, 
                   std::unique_ptr<CostFunctionInterface>& cost,
                   std::unique_ptr<ConstraintsInterface>& constraints) 
  : cost_(std::move(cost)),
    constraints_(std::move(constraints)),
    joint_constraints_(robot),
    riccati_matrix_factorizer_(robot),
    riccati_matrix_inverter_(robot),
    has_floating_base_(robot.has_floating_base()),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dim_passive_(robot.dim_passive()),
    max_dimf_(robot.max_dimf()),
    max_dimc_(robot.dim_passive()+robot.max_dimf()),
    dimf_(0),
    dimc_(0),
    lq_(Eigen::VectorXd::Zero(robot.dimv())),
    lv_(Eigen::VectorXd::Zero(robot.dimv())),
    la_(Eigen::VectorXd::Zero(robot.dimv())),
    lf_(Eigen::VectorXd::Zero(robot.max_dimf())),
    lu_(Eigen::VectorXd::Zero(robot.dimv())),
    lu_condensed_(Eigen::VectorXd::Zero(robot.dimv())),
    ka_(Eigen::VectorXd::Zero(robot.dimv())),
    kf_(Eigen::VectorXd::Zero(robot.max_dimf())),
    kmu_(Eigen::VectorXd::Zero(robot.max_dimf()+robot.dim_passive())),
    da_(Eigen::VectorXd::Zero(robot.dimv())),
    df_(Eigen::VectorXd::Zero(robot.max_dimf())),
    dmu_(Eigen::VectorXd::Zero(robot.max_dimf()+robot.dim_passive())),
    q_res_(Eigen::VectorXd::Zero(robot.dimv())),
    v_res_(Eigen::VectorXd::Zero(robot.dimv())),
    u_res_(Eigen::VectorXd::Zero(robot.dimv())),
    du_(Eigen::VectorXd::Zero(robot.dimv())),
    C_res_(Eigen::VectorXd::Zero(robot.max_dimf()+robot.dim_passive())),
    luu_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    du_dq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    du_dv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    du_da_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    du_df_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    Qqq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qqv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qqa_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qqf_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    Qvq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qvv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qva_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qvf_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    Qaa_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qaf_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    Qff_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.max_dimf())),
    Cq_(Eigen::MatrixXd::Zero(robot.max_dimf()+robot.dim_passive(), 
                              robot.dimv())),
    Cv_(Eigen::MatrixXd::Zero(robot.max_dimf()+robot.dim_passive(), 
                              robot.dimv())),
    Ca_(Eigen::MatrixXd::Zero(robot.max_dimf()+robot.dim_passive(), 
                              robot.dimv())),
    Cf_(Eigen::MatrixXd::Zero(robot.dim_passive(), robot.max_dimf())),
    Kaq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Kav_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Kfq_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
    Kfv_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
    Kmuq_(Eigen::MatrixXd::Zero(robot.max_dimf()+robot.dim_passive(),
                                robot.dimv())),
    Kmuv_(Eigen::MatrixXd::Zero(robot.max_dimf()+robot.dim_passive(), 
                                robot.dimv())),
    // The following variables are only needed for line search
    q_tmp_(Eigen::VectorXd::Zero(robot.dimq())), 
    v_tmp_(Eigen::VectorXd::Zero(robot.dimv())), 
    a_tmp_(Eigen::VectorXd::Zero(robot.dimv())), 
    f_tmp_(Eigen::VectorXd::Zero(robot.max_dimf())), 
    u_tmp_(Eigen::VectorXd::Zero(robot.dimv())), 
    u_res_tmp_(Eigen::VectorXd::Zero(robot.dimv())) {
}


SplitOCP::SplitOCP() 
  : cost_(),
    constraints_(),
    joint_constraints_(),
    riccati_matrix_factorizer_(),
    riccati_matrix_inverter_(),
    has_floating_base_(false),
    dimq_(0),
    dimv_(0),
    dimf_(0),
    max_dimf_(0),
    dimc_(0),
    dim_passive_(0),
    lq_(),
    lv_(),
    la_(),
    lf_(),
    lu_(),
    lu_condensed_(),
    ka_(),
    kf_(),
    kmu_(),
    da_(),
    df_(),
    dmu_(),
    q_res_(),
    v_res_(),
    u_res_(),
    du_(),
    C_res_(),
    luu_(),
    du_dq_(),
    du_dv_(),
    du_da_(),
    du_df_(),
    Qqq_(),
    Qqv_(),
    Qqa_(),
    Qqf_(),
    Qvq_(),
    Qvv_(),
    Qva_(),
    Qvf_(),
    Qaa_(),
    Qaf_(),
    Qff_(),
    Cq_(),
    Cv_(),
    Ca_(),
    Cf_(),
    Kaq_(),
    Kav_(),
    Kfq_(),
    Kfv_(),
    Kmuq_(),
    Kmuv_(),
    q_tmp_(), 
    v_tmp_(), 
    a_tmp_(), 
    f_tmp_(), 
    u_tmp_(), 
    u_res_tmp_() {
}


SplitOCP::~SplitOCP() {
}


bool SplitOCP::isFeasible(const Robot& robot, const Eigen::VectorXd& q, 
                          const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                          const Eigen::VectorXd& u) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(u.size() == dimv_);
  return joint_constraints_.isFeasible(q, v, a, u);
}


void SplitOCP::initConstraints(const Robot& robot, const int time_step, 
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


void SplitOCP::linearizeOCP(Robot& robot, const double t, const double dtau, 
                            const Eigen::VectorXd& lmd, 
                            const Eigen::VectorXd& gmm, 
                            const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                            const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                            const Eigen::VectorXd& f, const Eigen::VectorXd& mu,
                            const Eigen::VectorXd& lmd_next, 
                            const Eigen::VectorXd& gmm_next, 
                            const Eigen::VectorXd& q_next, 
                            const Eigen::VectorXd& v_next) {
  assert(dtau > 0);
  assert(lmd.size() == dimv_);
  assert(gmm.size() == dimv_);
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(u.size() == dimv_);
  assert(f.size() == max_dimf_);
  assert(mu.size() == max_dimc_);
  assert(lmd_next.size() == dimv_);
  assert(gmm_next.size() == dimv_);
  assert(q_next.size() == dimq_);
  assert(v_next.size() == dimv_);
  dimf_ = robot.dimf();
  dimc_ = robot.dim_passive() + robot.dimf();
  if (dimf_ > 0) {
    robot.updateKinematics(q, v, a);
  }
  ocplinearizer::linearizeStageCost(robot, cost_, t, dtau, q, v, a, u, f, 
                                    lq_, lv_, la_, lu_, lf_);
  ocplinearizer::linearizeDynamics(robot, dtau, q, v, a, u, f, q_next, v_next, 
                                   q_res_, v_res_, u_res_, du_dq_, du_dv_, 
                                   du_da_, du_df_);
  ocplinearizer::linearizeConstraints(robot, dtau, q, v, a, u, u_res_,  
                                      du_dq_, du_dv_, du_da_, du_df_, 
                                      C_res_, Cq_, Cv_, Ca_, Cf_);
  // Condense the control input torques and the Lagrange multiplier with 
  // respect to inverse dynamics.
  joint_constraints_.augmentDualResidual(dtau, lu_);
  cost_->luu(robot, t, dtau, u, luu_);
  joint_constraints_.condenseSlackAndDual(dtau, u, luu_, lu_);
  lu_condensed_ = lu_ + luu_ * u_res_;
  // Augment the condensed Newton residual of the contorl input torques. 
  lq_.noalias() += du_dq_.transpose() * lu_condensed_;
  lv_.noalias() += du_dv_.transpose() * lu_condensed_;
  la_.noalias() += du_da_.transpose() * lu_condensed_;
  lf_.head(dimf_).noalias() 
                += du_df_.leftCols(dimf_).transpose() * lu_condensed_;
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
  lf_.head(dimf_).noalias() += Cf_.leftCols(dimf_).transpose() 
      * mu.head(dim_passive_);
  // Augment the condensed Hessian of the contorl input torques. 
  Qqq_ = du_dq_.transpose() * luu_ * du_dq_;
  Qqv_ = du_dq_.transpose() * luu_ * du_dv_;
  Qqa_ = du_dq_.transpose() * luu_ * du_da_;
  Qvv_ = du_dv_.transpose() * luu_ * du_dv_;
  Qva_ = du_dv_.transpose() * luu_ * du_da_;
  Qaa_ = du_da_.transpose() * luu_ * du_da_;
  Qvq_ = Qqv_.transpose();
  Qqf_.leftCols(dimf_) = du_dq_.transpose() * luu_ * du_df_.leftCols(dimf_);
  Qvf_.leftCols(dimf_) = du_dv_.transpose() * luu_ * du_df_.leftCols(dimf_);
  Qaf_.leftCols(dimf_) = du_da_.transpose() * luu_ * du_df_.leftCols(dimf_);
  Qff_.topLeftCorner(dimf_, dimf_) 
      = du_df_.leftCols(dimf_).transpose() * luu_ * du_df_.leftCols(dimf_);
  // Condense the slack and dual variables of the inequality constraints on 
  // the configuration, velocity, and acceleration.
  joint_constraints_.condenseSlackAndDual(dtau, q, v, a, Qqq_, Qvv_, Qaa_, 
                                          lq_, lv_, la_);
  // Augment the cost function Hessian. 
  cost_->augment_lqq(robot, t, dtau, q, v, a, Qqq_);
  cost_->augment_lvv(robot, t, dtau, q, v, a, Qvv_);
  cost_->augment_laa(robot, t, dtau, q, v, a, Qaa_);
  cost_->augment_lff(robot, t, dtau, f, Qff_);
  if (robot.has_floating_base()) {
    riccati_matrix_factorizer_.setIntegrationSensitivities(robot, dtau, q, v);
  }
  if (dimf_ > 0) {
    riccati_matrix_inverter_.setContactStatus(robot);
    riccati_matrix_inverter_.precompute(Qaf_, Qff_);
  }
}


void SplitOCP::backwardRiccatiRecursion(const double dtau, 
                                        const Eigen::MatrixXd& Pqq_next, 
                                        const Eigen::MatrixXd& Pqv_next, 
                                        const Eigen::MatrixXd& Pvq_next, 
                                        const Eigen::MatrixXd& Pvv_next, 
                                        const Eigen::VectorXd& sq_next, 
                                        const Eigen::VectorXd& sv_next, 
                                        Eigen::MatrixXd& Pqq, 
                                        Eigen::MatrixXd& Pqv, 
                                        Eigen::MatrixXd& Pvq, 
                                        Eigen::MatrixXd& Pvv, 
                                        Eigen::VectorXd& sq, 
                                        Eigen::VectorXd& sv) {
  assert(dtau > 0);
  assert(Pqq_next.rows() == dimv_);
  assert(Pqq_next.cols() == dimv_);
  assert(Pqv_next.rows() == dimv_);
  assert(Pqv_next.cols() == dimv_);
  assert(Pvq_next.rows() == dimv_);
  assert(Pvq_next.cols() == dimv_);
  assert(Pvv_next.rows() == dimv_);
  assert(Pvv_next.cols() == dimv_);
  assert(sq_next.size() == dimv_);
  assert(sv_next.size() == dimv_);
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
  // Qqq_, Qqv_, Qvq_, Qvv_: representing Riccati factorization F
  // Qqa_, Qqf_, Qva_, Qvf_ : representing Riccati factorization H
  // Qaa_, Qaf_, Qff_ : representing Riccati factorization G
  riccati_matrix_factorizer_.factorize(dtau, Pqq_next, Pqv_next, Pvq_next, 
                                       Pvv_next, Qqq_, Qqv_, Qvq_, Qvv_);
  riccati_matrix_factorizer_.factorize(dtau, Pqv_next, Pvv_next, Qqa_, Qva_);
  riccati_matrix_factorizer_.factorize(dtau, Pvv_next, Qaa_);
  la_.noalias() += dtau * Pvq_next * q_res_;
  la_.noalias() += dtau * Pvv_next * v_res_;
  la_.noalias() -= dtau * sv_next;
  // Computes the state feedback gain and feedforward terms
  if (has_floating_base_) {
    if (dimf_ == 0) {
      riccati_matrix_inverter_.invert(Qqa_, Qva_, Qaa_, Cq_, Cv_, Ca_, la_, 
                                      C_res_, Kaq_, Kav_, Kmuq_, Kmuv_, ka_, 
                                      kmu_);
    }
    else if (dimf_ > 0) {
      riccati_matrix_inverter_.invert(Qqa_, Qva_, Qaa_, Qqf_, Qvf_, Cq_, Cv_, 
                                      Ca_, Cf_, la_, lf_, C_res_, Kaq_, Kav_, 
                                      Kfq_, Kfv_, Kmuq_, Kmuv_, ka_, kf_, kmu_);
    }
  } 
  else {
    if (dimf_ == 0) {
      riccati_matrix_inverter_.invert(Qqa_, Qva_, Qaa_, la_, Kaq_, Kav_, ka_);
    }
    else if (dimf_ > 0) {
      riccati_matrix_inverter_.invert(Qqa_, Qva_, Qaa_, Qqf_, Qvf_, Cq_, Cv_,  
                                      Ca_, la_, lf_, C_res_, Kaq_, Kav_, Kfq_,  
                                      Kfv_, Kmuq_, Kmuv_, ka_, kf_, kmu_);
    }
  }
  // Computes the Riccati factorization matrices
  Pqq = Qqq_;
  Pqq.noalias() += Kaq_.transpose() * Qqa_.transpose();
  Pqv = Qqv_;
  Pqv.noalias() += Kaq_.transpose() * Qva_.transpose();
  Pvq = Qvq_;
  Pvq.noalias() += Kav_.transpose() * Qqa_.transpose();
  Pvv = Qvv_;
  Pvv.noalias() += Kav_.transpose() * Qva_.transpose();
  Pqq.noalias() += Kfq_.topRows(dimf_).transpose() 
                  * Qqf_.leftCols(dimf_).transpose();
  Pqv.noalias() += Kfq_.topRows(dimf_).transpose() 
                  * Qvf_.leftCols(dimf_).transpose();
  Pvq.noalias() += Kfv_.topRows(dimf_).transpose() 
                  * Qqf_.leftCols(dimf_).transpose();
  Pvv.noalias() += Kfv_.topRows(dimf_).transpose() 
                  * Qvf_.leftCols(dimf_).transpose();
  Pqq.noalias() += Kmuq_.topRows(dimc_).transpose() * Cq_.topRows(dimc_);
  Pqv.noalias() += Kmuq_.topRows(dimc_).transpose() * Cv_.topRows(dimc_);
  Pvq.noalias() += Kmuv_.topRows(dimc_).transpose() * Cq_.topRows(dimc_);
  Pvv.noalias() += Kmuv_.topRows(dimc_).transpose() * Cv_.topRows(dimc_);
  // Computes the Riccati factorization vectors
  sq = sq_next - lq_;
  sq.noalias() -= Pqq_next * q_res_;
  sq.noalias() -= Pqv_next * v_res_;
  sq.noalias() -= Qqa_ * ka_;
  sv = dtau * sq_next + sv_next - lv_;
  sv.noalias() -= dtau * Pqq_next * q_res_;
  sv.noalias() -= Pvq_next * q_res_;
  sv.noalias() -= dtau * Pqv_next * v_res_;
  sv.noalias() -= Pvv_next * v_res_;
  sv.noalias() -= Qva_ * ka_;
  sq.noalias() -= Qqf_.leftCols(dimf_) * kf_.head(dimf_);
  sv.noalias() -= Qvf_.leftCols(dimf_) * kf_.head(dimf_);
  sq.noalias() -= Cq_.topRows(dimc_).transpose() * kmu_.head(dimc_);
  sv.noalias() -= Cv_.topRows(dimc_).transpose() * kmu_.head(dimc_);
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
  df_.head(dimf_) = kf_.head(dimf_) + Kfq_.topRows(dimf_) * dq 
                                    + Kfv_.topRows(dimf_) * dv;
  dmu_.head(dimc_) = kmu_.head(dimc_) + Kmuq_.topRows(dimc_) * dq 
                                      + Kmuv_.topRows(dimc_) * dv;
  du_ = u_res_;
  du_.noalias() += du_dq_ * dq;
  du_.noalias() += du_dv_ * dv;
  du_.noalias() += du_da_ * da_;
  du_.noalias() += du_df_.leftCols(dimf_) * df_.head(dimf_);
  joint_constraints_.computeSlackAndDualDirection(dtau, dq, dv, da_, du_);
}

 
double SplitOCP::maxPrimalStepSize() {
  return joint_constraints_.maxSlackStepSize();
}


double SplitOCP::maxDualStepSize() {
  return joint_constraints_.maxDualStepSize();
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
  f_tmp_.head(dimf_) = f.head(dimf_) + step_size * df_.head(dimf_);
  robot.setContactForces(f_tmp_);
  double cost = 0;
  cost += cost_->l(robot, t, dtau, q_tmp_, v_tmp_, a_tmp_, u_tmp_, f_tmp_);
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
  joint_constraints_.updateDual(step_size);
}


void SplitOCP::updatePrimal(Robot& robot, const double step_size, 
                            const double dtau, const Eigen::VectorXd& dq, 
                            const Eigen::VectorXd& dv, 
                            const Eigen::MatrixXd& Pqq, 
                            const Eigen::MatrixXd& Pqv, 
                            const Eigen::MatrixXd& Pvq, 
                            const Eigen::MatrixXd& Pvv, 
                            const Eigen::VectorXd& sq, 
                            const Eigen::VectorXd& sv, Eigen::VectorXd& q, 
                            Eigen::VectorXd& v, Eigen::VectorXd& a, 
                            Eigen::VectorXd& u, Eigen::VectorXd& beta, 
                            Eigen::VectorXd& f, Eigen::VectorXd& mu, 
                            Eigen::VectorXd& lmd, Eigen::VectorXd& gmm) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  assert(dq.size() == dimv_);
  assert(dv.size() == dimv_);
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
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(u.size() == dimv_);
  assert(beta.size() == dimv_);
  assert(f.size() == max_dimf_);
  assert(mu.size() == max_dimc_);
  assert(gmm.size() == dimv_);
  robot.integrateConfiguration(dq, step_size, q);
  v.noalias() += step_size * dv;
  a.noalias() += step_size * da_;
  u.noalias() += step_size * du_;
  f.head(dimf_).noalias() += step_size * df_.head(dimf_);
  mu.head(dimc_).noalias() += step_size * dmu_.head(dimc_);
  lu_.noalias() -= dtau * beta;
  beta.noalias() += step_size * lu_ / dtau;
  beta.noalias() += step_size * luu_ * du_ / dtau;
  lmd.noalias() += step_size * (Pqq * dq + Pqv * dv - sq);
  gmm.noalias() += step_size * (Pvq * dq + Pvv * dv - sv);
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
                                     const double dtau, 
                                     const Eigen::VectorXd& lmd, 
                                     const Eigen::VectorXd& gmm, 
                                     const Eigen::VectorXd& q, 
                                     const Eigen::VectorXd& v, 
                                     const Eigen::VectorXd& a, 
                                     const Eigen::VectorXd& u, 
                                     const Eigen::VectorXd& beta, 
                                     const Eigen::VectorXd& f, 
                                     const Eigen::VectorXd& mu, 
                                     const Eigen::VectorXd& lmd_next, 
                                     const Eigen::VectorXd& gmm_next, 
                                     const Eigen::VectorXd& q_next,
                                     const Eigen::VectorXd& v_next) {
  assert(dtau > 0);
  assert(lmd.size() == dimv_);
  assert(gmm.size() == dimv_);
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(u.size() == dimv_);
  assert(beta.size() == dimv_);
  assert(f.size() == max_dimf_);
  assert(mu.size() == max_dimc_);
  assert(lmd_next.size() == dimv_);
  assert(gmm_next.size() == dimv_);
  assert(q_next.size() == dimq_);
  assert(v_next.size() == dimv_);
  const int dimf = robot.dimf();
  const int dimc = robot.dimf() + robot.dim_passive();
  if (dimf > 0) {
    robot.updateKinematics(q, v, a);
  }
  ocplinearizer::linearizeStageCost(robot, cost_, t, dtau, q, v, a, u, f, 
                                    lq_, lv_, la_, lu_, lf_);
  ocplinearizer::linearizeDynamics(robot, dtau, q, v, a, u, f, q_next, v_next, 
                                   q_res_, v_res_, u_res_, du_dq_, du_dv_, 
                                   du_da_, du_df_);
  ocplinearizer::linearizeConstraints(robot, dtau, q, v, a, u, u_res_,  
                                      du_dq_, du_dv_, du_da_, du_df_, 
                                      C_res_, Cq_, Cv_, Ca_, Cf_);
  // Augment the dynamics constraints. 
  lq_.noalias() += lmd_next - lmd;
  lv_.noalias() += dtau * lmd_next + gmm_next - gmm;
  la_.noalias() += dtau * gmm_next;
  // Augment the partial derivatives of the inverse dynamics constraint. 
  lq_.noalias() += dtau * du_dq_.transpose() * beta;
  lv_.noalias() += dtau * du_dv_.transpose() * beta;
  la_.noalias() += dtau * du_da_.transpose() * beta;
  lu_.noalias() -= dtau * beta;
  lf_.head(dimf).noalias() += dtau * du_df_.leftCols(dimf).transpose() * beta;
  // Augmnet the partial derivatives of the inequality constriants.
  joint_constraints_.augmentDualResidual(dtau, lu_);
  joint_constraints_.augmentDualResidual(dtau, lq_, lv_, la_);
  // Augment the equality constraints 
  lq_.noalias() += Cq_.topRows(dimc).transpose() * mu.head(dimc);
  lv_.noalias() += Cv_.topRows(dimc).transpose() * mu.head(dimc);
  la_.noalias() += Ca_.topRows(dimc).transpose() * mu.head(dimc);
  lf_.head(dimf).noalias() += Cf_.leftCols(dimf).transpose() 
                                * mu.head(robot.dim_passive());
  double error = 0;
  error += q_res_.squaredNorm();
  error += v_res_.squaredNorm();
  error += u_res_.squaredNorm();
  error += lq_.squaredNorm();
  error += lv_.squaredNorm();
  error += la_.squaredNorm();
  error += lu_.squaredNorm();
  error += lf_.head(dimf).squaredNorm();
  error += joint_constraints_.residualSquaredNrom(dtau, q, v, a, u);
  error += C_res_.head(dimc).squaredNorm();
  return error;
}

} // namespace idocp