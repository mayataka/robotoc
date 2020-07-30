#include "ocp/split_ocp.hpp"

#include <assert.h>
#include "Eigen/LU"


namespace idocp {

SplitOCP::SplitOCP(const Robot& robot, const CostFunctionInterface* cost, 
                   const ConstraintsInterface* constraints) 
  : cost_(const_cast<CostFunctionInterface*>(cost)),
    constraints_(const_cast<ConstraintsInterface*>(constraints)),
    joint_constraints_(robot),
    riccati_matrix_factorizer_(robot),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimf_(0),
    dim_passive_(robot.dim_passive()),
    f_(Eigen::VectorXd::Zero(robot.max_dimf())),
    mu_(Eigen::VectorXd::Zero(robot.max_dimf()+robot.dim_passive())),
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
    a_res_(Eigen::VectorXd::Zero(robot.dimv())),
    f_res_(Eigen::VectorXd::Zero(robot.max_dimf())),
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
    Cf_(Eigen::MatrixXd::Zero(robot.max_dimf()+robot.dim_passive(), 
                              robot.max_dimf())),
    Qff_inv_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.max_dimf())),
    Saa_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Saa_inv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Saf_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    D_hat_(Eigen::MatrixXd::Zero(robot.max_dimf()+robot.dim_passive(), 
                                 robot.max_dimf()+robot.dim_passive())),
    D_hat_inv_(Eigen::MatrixXd::Zero(robot.max_dimf()+robot.dim_passive(), 
                                     robot.max_dimf()+robot.dim_passive())),
    L_U_(Eigen::MatrixXd::Zero(robot.dimv(), 
                               robot.max_dimf()+robot.dim_passive())),
    L_L_(Eigen::MatrixXd::Zero(robot.max_dimf(), 
                               robot.max_dimf()+robot.dim_passive())),
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


void SplitOCP::initConstraints(Robot& robot, const int time_step,
                               const double dtau, const Eigen::VectorXd& q, 
                               const Eigen::VectorXd& v, 
                               const Eigen::VectorXd& a, 
                               const Eigen::VectorXd& u) {
  assert(time_step >= 0);
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(a.size() == robot.dimv());
  assert(u.size() == robot.dimv());
  joint_constraints_.setTimeStep(time_step);
  joint_constraints_.setSlackAndDual(robot, dtau, q, v, a, u);
}


void SplitOCP::linearizeOCP(Robot& robot, const double t, const double dtau, 
                            const Eigen::VectorXd& lmd, 
                            const Eigen::VectorXd& gmm, 
                            const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                            const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
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
  assert(lmd_next.size() == robot.dimv());
  assert(gmm_next.size() == robot.dimv());
  assert(q_next.size() == robot.dimq());
  assert(v_next.size() == robot.dimv());
  // Residual of the state equation
  q_res_ = q + dtau * v - q_next;
  v_res_ = v + dtau * a - v_next;
  // First, we condense the the control input torques and the Lagrange 
  // multiplier with respect to inverse dynamics.
  // Partial derivatives of the cost function with respect to the control input
  // torques.
  cost_->lu(robot, t, dtau, u, lu_);
  // Augment the partial derivatives of the inequality constraint.
  joint_constraints_.augmentDualResidual(robot, dtau, lu_);
  // Hessian of the cost function.
  cost_->luu(robot, t, dtau, u, luu_);
  // Modify the Hessian and residual by condensing the slack and dual variables 
  // of the inequality constraints on the control input.
  joint_constraints_.condenseSlackAndDual(robot, dtau, u, luu_, lu_);
  // Get the present dimension of the contacts
  dimf_ = robot.dimf();
  if (dimf_ > 0) {
    robot.setContactForces(f_);
  }
  // Residual of the inverse dynamics constraint.
  robot.RNEA(q, v, a, u_res_);
  u_res_.noalias() -= u;
  // Condensed Newton residual with respect to the control input torques.
  lu_condensed_ = lu_ + luu_ * u_res_;
  // Partial derivatives of the cost function with respect to the configuration,
  // velocity, and acceleration.
  cost_->lq(robot, t, dtau, q, v, a, lq_);
  cost_->lv(robot, t, dtau, q, v, a, lv_);
  cost_->la(robot, t, dtau, q, v, a, la_);
  if (dimf_ > 0) {
    cost_->lf(robot, t, dtau, f_, lf_);
  }
  else {
  }
  // Augmnet the partial derivatives of the state equation.
  lq_.noalias() += lmd_next - lmd;
  lv_.noalias() += dtau * lmd_next + gmm_next - gmm;
  la_.noalias() += dtau * gmm_next;
  // Augmnet the partial derivatives of the inequality constriants.
  joint_constraints_.augmentDualResidual(robot, dtau, lq_, lv_, la_);
  // Augment the condensed Newton residual of the contorl input torques. 
  robot.RNEADerivatives(q, v, a, du_dq_, du_dv_, du_da_);
  lq_.noalias() += du_dq_.transpose() * lu_condensed_;
  lv_.noalias() += du_dv_.transpose() * lu_condensed_;
  la_.noalias() += du_da_.transpose() * lu_condensed_;
  if (dimf_ > 0) {
    // Condensing the input torque in the contact forces
    robot.updateKinematics(q, v, a);
    robot.dRNEAPartialdFext(du_df_);
    lf_.noalias() += du_df_.leftCols(dimf_).transpose() * lu_condensed_;
    // Computes the contact constraints.
    robot.computeBaumgarteResidual(C_res_);
    robot.computeBaumgarteDerivatives(Cq_, Cv_, Ca_);
    // Augment the equality constraints 
    lq_.noalias() += Cq_.topRows(dimf_).transpose() * mu_.head(dimf_);
    lv_.noalias() += Cv_.topRows(dimf_).transpose() * mu_.head(dimf_);
    la_.noalias() += Ca_.topRows(dimf_).transpose() * mu_.head(dimf_);
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
  // Modify the Hessian and residual by condensing the slack and dual variables 
  // of the inequality constraints on the configuration, velocity, and 
  // acceleration.
  joint_constraints_.condenseSlackAndDual(robot, dtau, q, v, a, Qqq_, Qvv_, 
                                          Qaa_, lq_, lv_, la_);
  // Augment the cost function Hessian. 
  cost_->lqq(robot, t, dtau, q, v, a, Qqq_);
  cost_->lvv(robot, t, dtau, q, v, a, Qvv_);
  cost_->laa(robot, t, dtau, q, v, a, Qaa_);
  if (dimf_ > 0) {
    cost_->lff(robot, t, dtau, f_, Qff_);
    // Precompute some parts of the Riccati factorization.
    Qff_inv_.topLeftCorner(dimf_, dimf_) 
        = Qff_.topLeftCorner(dimf_, dimf_)
              .llt().solve(Eigen::MatrixXd::Identity(dimf_, dimf_));
    Saa_ = - Qaf_.leftCols(dimf_) * Qff_inv_.topLeftCorner(dimf_, dimf_) 
                                  * Qaf_.leftCols(dimf_).transpose();
    Saf_.leftCols(dimf_) 
        = Qaf_.leftCols(dimf_) * Qff_inv_.topLeftCorner(dimf_, dimf_);
  }
  factorizer.computeIntegrationSensitivities(robot, dtau, q, v);
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
  assert(Qq.size() == robot.dimv());
  assert(Qv.size() == robot.dimv());
  cost_->phiq(robot, t, q, v, lq_);
  cost_->phiv(robot, t, q, v, lv_);
  Qq = - lq_ + lmd;
  Qv = - lv_ + gmm;
  cost_->phiqq(robot, t, q, v, Qqq);
  cost_->phivv(robot, t, q, v, Qvv);
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
  // Qqq_, Qqv_, Qvq_, Qvv_: representing Riccati factorizor F
  // Qqa_, Qqf_, Qva_, Qvf_ : representing Riccati factorizor H
  // Qaa_, Qaf_, Qff_ : representing Riccati factorizor G
  riccati_matrix_factorizer_.factorize(dtau, Pqq_next, Pqv_next, Pvq_next, 
                                       Pvv_next, Qqq_, Qqv_, Qvq_, Qvv_);
  riccati_matrix_factorizer_.factorize(dtau,  Pqv_next, Pvv_next, Qqa_, Qva_);
  riccati_matrix_factorizer_.factorize(dtau,  Pvv_next, Qaa_);
  if (dimf_ == 0) {
    // Computes inversion of the coefficient matrix of the decision variables
    Saa_inv_ = Qaa_.llt().solve(Eigen::MatrixXd::Identity(dimv_, dimv_));
    // Computes the state feedback gain
    Kaq_ = - Saa_inv_ * Qqa_.transpose();
    Kav_ = - Saa_inv_ * Qva_.transpose();
    // Computes the Riccati factorization matrices
    Pqq = Qqq_;
    Pqq.noalias() += Kaq_.transpose() * Qqa_.transpose();
    Pqv = Qqv_;
    Pqv.noalias() += Kaq_.transpose() * Qva_.transpose();
    Pvq = Qvq_;
    Pvq.noalias() += Kav_.transpose() * Qqa_.transpose();
    Pvv = Qvv_;
    Pvv.noalias() += Kav_.transpose() * Qva_.transpose();
    // Computes the feedforward terms
    la_.noalias() += dtau * Pvq_next * q_res_;
    la_.noalias() += dtau * Pvv_next * v_res_;
    la_.noalias() -= dtau * sv_next;
    ka_ = - Saa_inv_ * la_;
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
  }
  else if (dimf_ > 0) {
    // Computes inversion of the coefficient matrix of the decision variables
    Saa_.noalias() += Qaa_;
    Saa_inv_ = Saa_.llt().solve(Eigen::MatrixXd::Identity(dimv_, dimv_));
    Saf_.leftCols(dimf_) = Saa_inv_ * Saf_.leftCols(dimf_);
    D_hat_.topLeftCorner(dimf_, dimf_) 
        = - Ca_.topRows(dimf_) * Saa_inv_ * Ca_.topRows(dimf_).transpose();
    D_hat_inv_.topLeftCorner(dimf_, dimf_) 
        = D_hat_.topLeftCorner(dimf_, dimf_)
                .llt().solve(Eigen::MatrixXd::Identity(dimf_, dimf_));
    L_U_.leftCols(dimf_) 
        = Saa_inv_ * Ca_.topRows(dimf_).transpose() 
                   * D_hat_inv_.topLeftCorner(dimf_, dimf_);
    L_L_.topLeftCorner(dimf_, dimf_) 
        = - Saf_.leftCols(dimf_).transpose() 
            * Ca_.topRows(dimf_).transpose() 
            * D_hat_inv_.topLeftCorner(dimf_, dimf_);
    Saa_inv_.noalias() 
        += L_U_.leftCols(dimf_) * D_hat_.topLeftCorner(dimf_, dimf_)
                                * L_U_.leftCols(dimf_).transpose();
    Qff_inv_.noalias() 
        += Saf_.leftCols(dimf_).transpose() * Saa_ * Saf_.leftCols(dimf_);
    Qff_inv_.noalias() 
        += L_L_.topLeftCorner(dimf_, dimf_) 
            * D_hat_.topLeftCorner(dimf_, dimf_) 
            * L_L_.topLeftCorner(dimf_, dimf_).transpose();
    Saf_.leftCols(dimf_).noalias() 
        -= L_U_.leftCols(dimf_) * D_hat_.topLeftCorner(dimf_, dimf_)
                                * L_L_.topLeftCorner(dimf_, dimf_).transpose();
    // Computes the state feedback gain
    Kaq_ = - Saa_inv_ * Qqa_.transpose();
    Kaq_.noalias() += Saf_.leftCols(dimf_) * Qqf_.leftCols(dimf_).transpose();
    Kaq_.noalias() += L_U_.leftCols(dimf_) * Cq_.topRows(dimf_);
    Kav_ = - Saa_inv_ * Qva_.transpose();
    Kav_.noalias() += Saf_.leftCols(dimf_) * Qvf_.leftCols(dimf_).transpose();
    Kav_.noalias() += L_U_.leftCols(dimf_) * Cv_.topRows(dimf_);
    Kfq_.topRows(dimf_) = Saf_.leftCols(dimf_).transpose() * Qqa_.transpose();
    Kfq_.topRows(dimf_).noalias() 
        -= Qff_inv_.topLeftCorner(dimf_, dimf_) 
            * Qqf_.leftCols(dimf_).transpose();
    Kfq_.topRows(dimf_).noalias() 
        += L_L_.topLeftCorner(dimf_, dimf_) * Cq_.topRows(dimf_);
    Kfv_.topRows(dimf_) = Saf_.leftCols(dimf_).transpose() * Qva_.transpose();
    Kfv_.topRows(dimf_).noalias() 
        -= Qff_inv_.topLeftCorner(dimf_, dimf_) 
            * Qvf_.leftCols(dimf_).transpose();
    Kfv_.topRows(dimf_).noalias() 
        += L_L_.topLeftCorner(dimf_, dimf_) * Cv_.topRows(dimf_);
    Kmuq_.topRows(dimf_) = L_U_.leftCols(dimf_).transpose() * Qqa_.transpose();
    Kmuq_.topRows(dimf_).noalias() 
        += L_L_.topLeftCorner(dimf_, dimf_).transpose() 
            * Qqf_.leftCols(dimf_).transpose();
    Kmuq_.topRows(dimf_).noalias() 
        -= D_hat_inv_.topLeftCorner(dimf_, dimf_) * Cq_.topRows(dimf_);
    Kmuv_.topRows(dimf_) = L_U_.leftCols(dimf_).transpose() * Qva_.transpose();
    Kmuv_.topRows(dimf_).noalias() 
        += L_L_.topLeftCorner(dimf_, dimf_).transpose() 
            * Qvf_.leftCols(dimf_).transpose();
    Kmuv_.topRows(dimf_).noalias() 
        -= D_hat_inv_.topLeftCorner(dimf_, dimf_) * Cv_.topRows(dimf_);
    // Computes the Riccati factorization matrices
    Pqq = Qqq_;
    Pqq.noalias() += Kaq_.transpose() * Qqa_.transpose();
    Pqq.noalias() += Kfq_.topRows(dimf_).transpose() 
                     * Qqf_.leftCols(dimf_).transpose();
    Pqq.noalias() += Kmuq_.topRows(dimf_).transpose() * Cq_.topRows(dimf_);
    Pqv = Qqv_;
    Pqv.noalias() += Kaq_.transpose() * Qva_.transpose();
    Pqv.noalias() += Kfq_.topRows(dimf_).transpose() 
                     * Qvf_.leftCols(dimf_).transpose();
    Pqv.noalias() += Kmuq_.topRows(dimf_).transpose() * Cv_.topRows(dimf_);
    Pvq = Qvq_;
    Pvq.noalias() += Kav_.transpose() * Qqa_.transpose();
    Pvq.noalias() += Kfv_.topRows(dimf_).transpose() 
                     * Qqf_.leftCols(dimf_).transpose();
    Pvq.noalias() += Kmuv_.topRows(dimf_).transpose() * Cq_.topRows(dimf_);
    Pvv = Qvv_;
    Pvv.noalias() += Kav_.transpose() * Qva_.transpose();
    Pvv.noalias() += Kfv_.topRows(dimf_).transpose() 
                     * Qvf_.leftCols(dimf_).transpose();
    Pvv.noalias() += Kmuv_.topRows(dimf_).transpose() * Cv_.topRows(dimf_);
    // Computes the feedforward terms
    la_.noalias() += dtau * Pvq_next * q_res_;
    la_.noalias() += dtau * Pvv_next * v_res_;
    la_.noalias() -= dtau * sv_next;
    ka_ = - Saa_inv_ * la_;
    ka_.noalias() += Saf_.leftCols(dimf_) * lf_.head(dimf_);
    ka_.noalias() += L_U_.leftCols(dimf_) * C_res_.head(dimf_);
    kf_.head(dimf_) = Saf_.leftCols(dimf_).transpose() * la_;
    kf_.head(dimf_).noalias() 
        -= Qff_inv_.topLeftCorner(dimf_, dimf_) * lf_.head(dimf_);
    kf_.head(dimf_).noalias() 
        += L_L_.topLeftCorner(dimf_, dimf_) * C_res_.head(dimf_);
    kmu_.head(dimf_) = L_U_.leftCols(dimf_).transpose() * la_;
    kmu_.head(dimf_).noalias() 
        += L_L_.topLeftCorner(dimf_, dimf_).transpose() * lf_.head(dimf_);
    kmu_.head(dimf_).noalias() 
        -= D_hat_inv_.topLeftCorner(dimf_, dimf_) * C_res_.head(dimf_);
    // Computes the Riccati factorization vectors
    sq = sq_next - lq_;
    sq.noalias() -= Pqq_next * q_res_;
    sq.noalias() -= Pqv_next * v_res_;
    sq.noalias() -= Qqa_ * ka_;
    sq.noalias() -= Qqf_.leftCols(dimf_) * kf_.head(dimf_);
    sq.noalias() -= Cq_.topRows(dimf_).transpose() * kmu_.head(dimf_);
    sv = dtau * sq_next + sv_next - lv_;
    sv.noalias() -= dtau * Pqq_next * q_res_;
    sv.noalias() -= Pvq_next * q_res_;
    sv.noalias() -= dtau * Pqv_next * v_res_;
    sv.noalias() -= Pvv_next * v_res_;
    sv.noalias() -= Qva_ * ka_;
    sv.noalias() -= Qvf_.leftCols(dimf_) * kf_.head(dimf_);
    sv.noalias() -= Cv_.topRows(dimf_).transpose() * kmu_.head(dimf_);
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


void SplitOCP::computeCondensedDirection(Robot& robot, const double dtau, 
                                         const Eigen::VectorXd& dq, 
                                         const Eigen::VectorXd& dv) {
  assert(dtau > 0);
  assert(dq.size() == robot.dimv());
  assert(dv.size() == robot.dimv());
  if (dimf_ > 0) {
    df_.head(dimf_) = kf_.head(dimf_) + Kfq_.topRows(dimf_) * dq 
                                      + Kfv_.topRows(dimf_) * dv;
    dmu_.head(dimf_) = kmu_.head(dimf_) + Kmuq_.topRows(dimf_) * dq 
                                        + Kmuv_.topRows(dimf_) * dv;
  }
  du_ = u_res_;
  du_.noalias() += du_dq_ * dq;
  du_.noalias() += du_dv_ * dv;
  du_.noalias() += du_da_ * da_;
  if (dimf_ > 0) {
    du_.noalias() += du_df_.leftCols(dimf_) * df_.head(dimf_);
  }
  joint_constraints_.computeSlackAndDualDirection(robot, dtau, dq, dv, da_, du_);
}

 
double SplitOCP::maxPrimalStepSize() {
  return joint_constraints_.maxSlackStepSize();
}


double SplitOCP::maxDualStepSize() {
  return joint_constraints_.maxDualStepSize();
}


double SplitOCP::stageCostDerivativeDotDirection(Robot& robot, const double t, 
                                                 const double dtau, 
                                                 const Eigen::VectorXd& q, 
                                                 const Eigen::VectorXd& v, 
                                                 const Eigen::VectorXd& a, 
                                                 const Eigen::VectorXd& u, 
                                                 const Eigen::VectorXd& dq, 
                                                 const Eigen::VectorXd& dv) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(a.size() == robot.dimv());
  assert(u.size() == robot.dimv());
  assert(dq.size() == robot.dimv());
  assert(dv.size() == robot.dimv());
  cost_->lq(robot, t, dtau, q, v, a, lq_);
  cost_->lv(robot, t, dtau, q, v, a, lv_);
  cost_->la(robot, t, dtau, q, v, a, la_);
  cost_->lu(robot, t, dtau, u, lu_);
  double product = 0;
  product += lq_.dot(dq);
  product += lv_.dot(dv);
  product += la_.dot(da_);
  product += lu_.dot(du_);
  return product;
}


double SplitOCP::terminalCostDerivativeDotDirection(Robot& robot, 
                                                    const double t, 
                                                    const Eigen::VectorXd& q, 
                                                    const Eigen::VectorXd& v, 
                                                    const Eigen::VectorXd& dq,
                                                    const Eigen::VectorXd& dv) {
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(dq.size() == robot.dimv());
  assert(dv.size() == robot.dimv());
  cost_->phiq(robot, t, q, v, lq_);
  cost_->phiv(robot, t, q, v, lv_);
  double product = 0;
  product += lq_.dot(dq);
  product += lv_.dot(dv);
  return product;
}


std::pair<double, double> SplitOCP::stageCostAndConstraintsViolation(
    Robot& robot, const double t, const double dtau, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
    const Eigen::VectorXd& a, const Eigen::VectorXd& u) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(a.size() == robot.dimv());
  assert(u.size() == robot.dimv());
  double cost = 0;
  cost += cost_->l(robot, t, dtau, q, v, a, u);
  cost += joint_constraints_.costSlackBarrier();
  double constraints_violation = 0;
  constraints_violation += q_res_.lpNorm<1>();
  constraints_violation += v_res_.lpNorm<1>();
  constraints_violation += dtau * u_res_.lpNorm<1>();
  constraints_violation += joint_constraints_.residualL1Nrom(robot, dtau, q, v, 
                                                             a, u);
  return std::make_pair(cost, constraints_violation);
}


std::pair<double, double> SplitOCP::stageCostAndConstraintsViolation(
    Robot& robot, const double step_size, const double t, const double dtau, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
    const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
    const Eigen::VectorXd& q_next, const Eigen::VectorXd& v_next, 
    const Eigen::VectorXd& dq, const Eigen::VectorXd& dv, 
    const Eigen::VectorXd& dq_next, const Eigen::VectorXd& dv_next) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(a.size() == robot.dimv());
  assert(u.size() == robot.dimv());
  assert(q_next.size() == robot.dimq());
  assert(v_next.size() == robot.dimv());
  assert(dq.size() == robot.dimv());
  assert(dv.size() == robot.dimv());
  assert(dq_next.size() == robot.dimv());
  assert(dv_next.size() == robot.dimv());
  q_tmp_ = q + step_size * dq;
  v_tmp_ = v + step_size * dv;
  a_tmp_ = a + step_size * da_;
  u_tmp_ = u + step_size * du_;
  double cost = 0;
  cost += cost_->l(robot, t, dtau, q_tmp_, v_tmp_, a_tmp_, u_tmp_);
  cost += joint_constraints_.costSlackBarrier(step_size);
  q_res_ = q_tmp_ + dtau * v_tmp_ - q_next - step_size * dq_next;
  v_res_ = v_tmp_ + dtau * a_tmp_ - v_next - step_size * dv_next;
  robot.RNEA(q_tmp_, v_tmp_, a_tmp_, u_res_tmp_);
  u_res_tmp_.noalias() -= u_tmp_;
  double constraints_violation = 0;
  constraints_violation += q_res_.lpNorm<1>();
  constraints_violation += v_res_.lpNorm<1>();
  constraints_violation += joint_constraints_.residualL1Nrom(robot, dtau, 
                                                             q_tmp_, v_tmp_, 
                                                             a_tmp_, u_tmp_);
  constraints_violation += dtau * u_res_tmp_.lpNorm<1>();
  return std::make_pair(cost, constraints_violation);
}


double SplitOCP::terminalCost(Robot& robot, const double t, 
                              const Eigen::VectorXd& q, 
                              const Eigen::VectorXd& v) {
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  return cost_->phi(robot, t, q, v);
}


double SplitOCP::terminalCost(Robot& robot, const double step_size, 
                              const double t, const Eigen::VectorXd& q, 
                              const Eigen::VectorXd& v, 
                              const Eigen::VectorXd& dq, 
                              const Eigen::VectorXd& dv) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(dq.size() == robot.dimv());
  assert(dv.size() == robot.dimv());
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
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(a.size() == robot.dimv());
  assert(lmd.size() == robot.dimv());
  assert(gmm.size() == robot.dimv());
  q.noalias() += step_size * dq;
  v.noalias() += step_size * dv;
  a.noalias() += step_size * da_;
  u.noalias() += step_size * du_;
  if (dimf_ > 0) {
    f_.head(dimf_).noalias() += step_size * df_;
    mu_.head(dimf_).noalias() += step_size * dmu_;
  }
  // modify lu_ so that it includes beta.
  lu_.noalias() -= dtau * beta;
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
                            Eigen::VectorXd& gmm) const {
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


void SplitOCP::getStateFeedbackGain(Eigen::MatrixXd& Kq, 
                                    Eigen::MatrixXd& Kv) const {
  assert(Kq.cols() == dimv_);
  assert(Kq.rows() == dimv_);
  assert(Kv.cols() == dimv_);
  assert(Kv.rows() == dimv_);
  if (dimf_ == 0) {
    Kq = du_dq_ + du_da_ * Kaq_;
    Kv = du_dv_ + du_da_ * Kav_;
  }
  else {
    Kq = du_dq_ + du_da_ * Kaq_ + du_df_.leftCols(dimf_) * Kfq_.topRows(dimf_);
    Kv = du_dv_ + du_da_ * Kav_ + du_df_.leftCols(dimf_) * Kfv_.topRows(dimf_);
  }
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
  // Compute the partial derivatives of the Lagrangian with respect to the 
  // configuration, velocity, acceleration, and the control input torques.
  // Partial derivatives of the cost function.
  const int dimf = robot.dimf();
  cost_->lq(robot, t, dtau, q, v, a, lq_);
  cost_->lv(robot, t, dtau, q, v, a, lv_);
  cost_->la(robot, t, dtau, q, v, a, la_);
  cost_->lu(robot, t, dtau, u, lu_);
  if (dimf > 0) {
    cost_->lf(robot, t, dtau, f_, lf_);
  }
  // Augment the partial derivatives of the state equation.
  lq_.noalias() += lmd_next - lmd;
  lv_.noalias() += dtau * lmd_next + gmm_next - gmm;
  la_.noalias() += dtau * gmm_next;
  // Augment the partial derivatives of the inverse dynamics constraint.
  robot.RNEADerivatives(q, v, a, du_dq_, du_dv_, du_da_);
  lq_.noalias() += dtau * du_dq_.transpose() * beta;
  lv_.noalias() += dtau * du_dv_.transpose() * beta;
  la_.noalias() += dtau * du_da_.transpose() * beta;
  lu_.noalias() -= dtau * beta;
  if (dimf > 0) {
    robot.updateKinematics(q, v, a);
    robot.dRNEAPartialdFext(du_df_);
    lf_.noalias() += du_df_.leftCols(dimf_).transpose() * beta;
    // Computes the contact constraints.
    robot.computeBaumgarteResidual(C_res_);
    robot.computeBaumgarteDerivatives(Cq_, Cv_, Ca_);
    // Augment the equality constraints 
    lq_.noalias() += Cq_.topRows(dimf_).transpose() * mu_.head(dimf_);
    lv_.noalias() += Cv_.topRows(dimf_).transpose() * mu_.head(dimf_);
    la_.noalias() += Ca_.topRows(dimf_).transpose() * mu_.head(dimf_);
  }
  // Augment the partial derivatives of the inequality constraint.
  joint_constraints_.augmentDualResidual(robot, dtau, lq_, lv_, la_);
  joint_constraints_.augmentDualResidual(robot, dtau, lu_);
  // Compute the residual of the state eqation.
  q_res_ = q + dtau * v - q_next;
  v_res_ = v + dtau * a - v_next;
  // Compute the residual of the inverse dynamics constraint.
  if (dimf > 0) {
    robot.setContactForces(f_);
  }
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
  if (dimf > 0) {
    error += lf_.head(dimf).squaredNorm();
    error += C_res_.head(dimf).squaredNorm();
  }
  return error;
}


double SplitOCP::squaredKKTErrorNorm(Robot& robot, const double t, 
                                     const Eigen::VectorXd& lmd, 
                                     const Eigen::VectorXd& gmm, 
                                     const Eigen::VectorXd& q, 
                                     const Eigen::VectorXd& v) {
  assert(lmd.size() == robot.dimv());
  assert(gmm.size() == robot.dimv());
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  // Compute the partial derivatives of the Lagrangian with respect to the 
  // terminal configuration and velocity.
  cost_->phiq(robot, t, q, v, lq_);
  cost_->phiv(robot, t, q, v, lv_);
  lq_.noalias() -= lmd;
  lv_.noalias() -= gmm;
  double error = 0;
  error += lq_.squaredNorm();
  error += lv_.squaredNorm();
  return error;
}

} // namespace idocp