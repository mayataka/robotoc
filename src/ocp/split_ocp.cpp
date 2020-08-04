#include "idocp/ocp/split_ocp.hpp"

#include <assert.h>


namespace idocp {

SplitOCP::SplitOCP(const Robot& robot, const CostFunctionInterface* cost, 
                   const ConstraintsInterface* constraints) 
  : cost_(const_cast<CostFunctionInterface*>(cost)),
    constraints_(const_cast<ConstraintsInterface*>(constraints)),
    joint_constraints_(robot),
    riccati_matrix_factorizer_(robot),
    riccati_matrix_inverter_(robot),
    has_floating_base_(robot.has_floating_base()),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimf_(0),
    max_dimf_(robot.max_dimf()),
    dimc_(0),
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
    has_floating_base_(),
    dimq_(0),
    dimv_(0),
    dimf_(0),
    max_dimf_(0),
    dimc_(0),
    dim_passive_(0),
    f_(),
    mu_(),
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


void SplitOCP::setCostFunction(const CostFunctionInterface* cost) {
  cost_ = const_cast<CostFunctionInterface*>(cost);
}


void SplitOCP::setConstraints(const ConstraintsInterface* constraints) {
  constraints_ = const_cast<ConstraintsInterface*>(constraints);
}


void SplitOCP::set_f(const Eigen::VectorXd& f) {
  assert(f.size() == max_dimf_);
  f_ = f;
}


void SplitOCP::set_mu(const Eigen::VectorXd& mu) {
  assert(mu.size() == max_dimf_+dim_passive_);
  mu_ = mu;
}


bool SplitOCP::isFeasible(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                          const Eigen::VectorXd& a, const Eigen::VectorXd& u) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(u.size() == dimv_);
  return joint_constraints_.isFeasible(q, v, a, u);
}


void SplitOCP::initConstraints(const int time_step, const double dtau, 
                               const Eigen::VectorXd& q, 
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
  assert(lmd_next.size() == dimv_);
  assert(gmm_next.size() == dimv_);
  assert(q_next.size() == dimq_);
  assert(v_next.size() == dimv_);
  // Residual of the state equation
  robot.differenceConfiguration(q, q_next, q_res_);
  q_res_.noalias() += dtau * v;
  v_res_ = v + dtau * a - v_next;
  // First, we condense the the control input torques and the Lagrange 
  // multiplier with respect to inverse dynamics.
  // Partial derivatives of the cost function with respect to the input torques.
  cost_->lu(t, dtau, u, lu_);
  // Augment the partial derivatives of the inequality constraint.
  joint_constraints_.augmentDualResidual(dtau, lu_);
  // Hessian of the cost function.
  cost_->luu(t, dtau, u, luu_);
  // Modify the Hessian and residual by condensing the slack and dual variables 
  // of the inequality constraints on the control input.
  joint_constraints_.condenseSlackAndDual(dtau, u, luu_, lu_);
  // Get the present dimension of the contacts
  dimf_ = robot.dimf();
  dimc_ = robot.dimf() + robot.dim_passive();
  // Set the contact forces to the robot model and contact status to the cost.
  if (dimf_ > 0) {
    robot.setContactForces(f_);
    cost_->setContactStatus(robot);
  }
  // Residual of the inverse dynamics constraint.
  robot.RNEA(q, v, a, u_res_);
  u_res_.noalias() -= u;
  // Compute condensed Newton residual with respect to the input torques.
  lu_condensed_ = lu_ + luu_ * u_res_;
  // Partial derivatives of the cost function with respect to the configuration,
  // velocity, and acceleration.
  // If the robot has a floating base, set the configuration Jacobian into cost.
  if (has_floating_base_) {
    cost_->setConfigurationJacobian(robot, q);
  }
  cost_->lq(t, dtau, q, v, a, lq_);
  cost_->lv(t, dtau, q, v, a, lv_);
  cost_->la(t, dtau, q, v, a, la_);
  if (dimf_ > 0) {
    // Partial derivatives of the cost function with respect to the contact 
    // forces.
    cost_->lf(t, dtau, f_, lf_);
  }
  // Augmnet the partial derivatives of the state equation.
  lq_.noalias() += lmd_next - lmd;
  lv_.noalias() += dtau * lmd_next + gmm_next - gmm;
  la_.noalias() += dtau * gmm_next;
  // Augmnet the partial derivatives of the inequality constriants.
  joint_constraints_.augmentDualResidual(dtau, lq_, lv_, la_);
  robot.RNEADerivatives(q, v, a, du_dq_, du_dv_, du_da_);
  if (dimf_ > 0) {
    // Partial derivatives of the cost function with respect to the contact 
    // forces.
    robot.updateKinematics(q, v, a);
    robot.dRNEAPartialdFext(du_df_);
  }
  // Augment the condensed Newton residual of the contorl input torques. 
  lq_.noalias() += du_dq_.transpose() * lu_condensed_;
  lv_.noalias() += du_dv_.transpose() * lu_condensed_;
  la_.noalias() += du_da_.transpose() * lu_condensed_;
  if (dimf_ > 0) {
    // Condensing the input torque in the contact forces
    lf_.head(dimf_).noalias() += du_df_.leftCols(dimf_).transpose() * lu_condensed_;
  }
  if (has_floating_base_) {
    // The equality constraints of the floating base.
    Cq_.topRows(dim_passive_) = dtau * du_dq_.topRows(dim_passive_);
    Cv_.topRows(dim_passive_) = dtau * du_dv_.topRows(dim_passive_);
    Ca_.topRows(dim_passive_) = dtau * du_da_.topRows(dim_passive_);
    if (dimf_ > 0) {
      Cf_.leftCols(dimf_) = dtau * du_df_.topLeftCorner(dim_passive_, dimf_);
    }
    C_res_.head(dim_passive_) = dtau * u.head(dim_passive_);
    C_res_.head(dim_passive_).noalias() += dtau * u_res_.head(dim_passive_);
  }
  if (dimf_ > 0) {
    // The equality constraints of the contact constraints.
    robot.computeBaumgarteResidual(dim_passive_, dtau, C_res_);
    robot.computeBaumgarteDerivatives(dim_passive_, dtau, Cq_, Cv_, Ca_);
  }
  if (dimc_ > 0) {
    // Augment the equality constraints 
    lq_.noalias() += Cq_.topRows(dimc_).transpose() * mu_.head(dimc_);
    lv_.noalias() += Cv_.topRows(dimc_).transpose() * mu_.head(dimc_);
    la_.noalias() += Ca_.topRows(dimc_).transpose() * mu_.head(dimc_);
    if (dimf_ > 0 && dim_passive_ > 0) {
      lf_.head(dimf_).noalias() += Cf_.leftCols(dimf_).transpose() 
          * mu_.head(dim_passive_);
    }
  }
  if (has_floating_base_) {
    riccati_matrix_factorizer_.setIntegrationSensitivities(robot, dtau, q, v);
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
  joint_constraints_.condenseSlackAndDual(dtau, q, v, a, Qqq_, Qvv_, Qaa_, 
                                          lq_, lv_, la_);
  // Augment the cost function Hessian. 
  cost_->augment_lqq(t, dtau, q, v, a, Qqq_);
  cost_->augment_lvv(t, dtau, q, v, a, Qvv_);
  cost_->augment_laa(t, dtau, q, v, a, Qaa_);
  if (dimf_ > 0) {
    cost_->augment_lff(t, dtau, f_, Qff_);
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
  if (dimf_ > 0) {
    Pqq.noalias() += Kfq_.topRows(dimf_).transpose() 
                    * Qqf_.leftCols(dimf_).transpose();
    Pqv.noalias() += Kfq_.topRows(dimf_).transpose() 
                    * Qvf_.leftCols(dimf_).transpose();
    Pvq.noalias() += Kfv_.topRows(dimf_).transpose() 
                    * Qqf_.leftCols(dimf_).transpose();
    Pvv.noalias() += Kfv_.topRows(dimf_).transpose() 
                    * Qvf_.leftCols(dimf_).transpose();
  } 
  if (dimc_ > 0) {
    Pqq.noalias() += Kmuq_.topRows(dimc_).transpose() * Cq_.topRows(dimc_);
    Pqv.noalias() += Kmuq_.topRows(dimc_).transpose() * Cv_.topRows(dimc_);
    Pvq.noalias() += Kmuv_.topRows(dimc_).transpose() * Cq_.topRows(dimc_);
    Pvv.noalias() += Kmuv_.topRows(dimc_).transpose() * Cv_.topRows(dimc_);
  }
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
  if (dimf_ > 0) {
    sq.noalias() -= Qqf_.leftCols(dimf_) * kf_.head(dimf_);
    sv.noalias() -= Qvf_.leftCols(dimf_) * kf_.head(dimf_);
  }
  if (dimc_ > 0) {
    sq.noalias() -= Cq_.topRows(dimc_).transpose() * kmu_.head(dimc_);
    sv.noalias() -= Cv_.topRows(dimc_).transpose() * kmu_.head(dimc_);
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
    const Eigen::VectorXd& q_next, const Eigen::VectorXd& v_next, 
    const Eigen::VectorXd& dq, const Eigen::VectorXd& dv, 
    const Eigen::VectorXd& dq_next, const Eigen::VectorXd& dv_next) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(u.size() == dimv_);
  assert(q_next.size() == dimq_);
  assert(v_next.size() == dimv_);
  assert(dq.size() == dimv_);
  assert(dv.size() == dimv_);
  assert(dq_next.size() == dimv_);
  assert(dv_next.size() == dimv_);
  if (has_floating_base_) {
    q_tmp_ = q;
    robot.integrateConfiguration(dq, step_size, q_tmp_);
    cost_->setConfigurationJacobian(robot, q_tmp_);
  }
  else {
    q_tmp_ = q + step_size * dq;
  }
  v_tmp_ = v + step_size * dv;
  a_tmp_ = a + step_size * da_;
  u_tmp_ = u + step_size * du_;
  if (dimf_ > 0) {
    f_tmp_.head(dimf_) = f_.head(dimf_) + step_size * df_.head(dimf_);
    robot.setContactForces(f_tmp_);
  }
  double cost = 0;
  cost += cost_->l(t, dtau, q_tmp_, v_tmp_, a_tmp_, u_tmp_, f_tmp_);
  cost += joint_constraints_.costSlackBarrier(step_size);
  robot.differenceConfiguration(q_tmp_, q_next, q_res_);
  q_res_.noalias() += dtau * v_tmp_ - step_size * dq_next;
  v_res_ = v_tmp_ + dtau * a_tmp_ - v_next - step_size * dv_next;
  robot.RNEA(q_tmp_, v_tmp_, a_tmp_, u_res_tmp_);
  u_res_tmp_.noalias() -= u_tmp_;
  double constraints_violation = 0;
  constraints_violation += q_res_.lpNorm<1>();
  constraints_violation += v_res_.lpNorm<1>();
  constraints_violation += joint_constraints_.residualL1Nrom(dtau, q_tmp_, 
                                                             v_tmp_, a_tmp_, 
                                                             u_tmp_);
  constraints_violation += dtau * u_res_tmp_.lpNorm<1>();
  if (dimc_ > 0) {
    C_res_.head(dim_passive_) = dtau * u_tmp_.head(dim_passive_);
    robot.updateKinematics(q_tmp_, v_tmp_, a_tmp_);
    robot.computeBaumgarteResidual(dim_passive_, dtau, C_res_);
    constraints_violation += C_res_.head(dimc_).lpNorm<1>();
  }
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
  assert(lmd.size() == dimv_);
  assert(gmm.size() == dimv_);
  robot.integrateConfiguration(dq, step_size, q);
  v.noalias() += step_size * dv;
  a.noalias() += step_size * da_;
  u.noalias() += step_size * du_;
  if (dimf_ > 0) {
    f_.head(dimf_).noalias() += step_size * df_.head(dimf_);
    mu_.head(dimc_).noalias() += step_size * dmu_.head(dimc_);
  }
  // modify lu_ so that it includes beta.
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
  assert(lmd.size() == dimv_);
  assert(gmm.size() == dimv_);
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(a.size() == dimv_);
  assert(u.size() == dimv_);
  assert(beta.size() == dimv_);
  assert(lmd_next.size() == dimv_);
  assert(gmm_next.size() == dimv_);
  assert(q_next.size() == dimq_);
  assert(v_next.size() == dimv_);
  // Compute the partial derivatives of the Lagrangian with respect to the 
  // configuration, velocity, acceleration, and the control input torques.
  // Partial derivatives of the cost function.
  const int dimf = robot.dimf();
  const int dimc = robot.dimf() + robot.dim_passive();
  // Residual of the state eqation
  robot.differenceConfiguration(q, q_next, q_res_);
  q_res_.noalias() += dtau * v;
  v_res_ = v + dtau * a - v_next;
  // Compute the residual of the inverse dynamics constraint.
  if (dimf > 0) {
    robot.setContactForces(f_);
    cost_->setContactStatus(robot);
  }
  robot.RNEA(q, v, a, u_res_);
  u_res_.noalias() -= u;
  u_res_.array() *= dtau;
  if (has_floating_base_) {
    cost_->setConfigurationJacobian(robot, q);
  }
  cost_->lq(t, dtau, q, v, a, lq_);
  cost_->lv(t, dtau, q, v, a, lv_);
  cost_->la(t, dtau, q, v, a, la_);
  cost_->lu(t, dtau, u, lu_);
  if (dimf > 0) {
    cost_->setContactStatus(robot);
    cost_->lf(t, dtau, f_, lf_);
  }
  // Augment the partial derivatives of the state equation.
  lq_.noalias() += lmd_next - lmd;
  lv_.noalias() += dtau * lmd_next + gmm_next - gmm;
  la_.noalias() += dtau * gmm_next;
  // Augment the partial derivatives of the inequality constraints.
  joint_constraints_.augmentDualResidual(dtau, lq_, lv_, la_);
  joint_constraints_.augmentDualResidual(dtau, lu_);
  // Augment the partial derivatives of the inverse dynamics constraint.
  robot.RNEADerivatives(q, v, a, du_dq_, du_dv_, du_da_);
  lq_.noalias() += dtau * du_dq_.transpose() * beta;
  lv_.noalias() += dtau * du_dv_.transpose() * beta;
  la_.noalias() += dtau * du_da_.transpose() * beta;
  lu_.noalias() -= dtau * beta;
  if (dimf > 0) {
    robot.updateKinematics(q, v, a);
    robot.dRNEAPartialdFext(du_df_);
    lf_.head(dimf).noalias() += dtau * du_df_.leftCols(dimf).transpose() * beta;
  }
  if (has_floating_base_) {
    // The equality constraints of the floating base.
    Cq_.topRows(dim_passive_) = dtau * du_dq_.topRows(dim_passive_);
    Cv_.topRows(dim_passive_) = dtau * du_dv_.topRows(dim_passive_);
    Ca_.topRows(dim_passive_) = dtau * du_da_.topRows(dim_passive_);
    if (dimf > 0) {
      Cf_.leftCols(dimf) = dtau * du_df_.topLeftCorner(dim_passive_, dimf);
    }
    C_res_.head(dim_passive_) = dtau * u.head(dim_passive_);
  }
  if (dimf > 0) {
    // Computes the contact constraints.
    robot.computeBaumgarteResidual(dim_passive_, dtau, C_res_);
    robot.computeBaumgarteDerivatives(dim_passive_, dtau, Cq_, Cv_, Ca_);
  }
  if (dimc > 0) {
    // Augment the equality constraints 
    lq_.noalias() += Cq_.topRows(dimc).transpose() * mu_.head(dimc);
    lv_.noalias() += Cv_.topRows(dimc).transpose() * mu_.head(dimc);
    la_.noalias() += Ca_.topRows(dimc).transpose() * mu_.head(dimc);
    if (dimf > 0) {
      lf_.head(dimf).noalias() 
          += Cf_.leftCols(dimf).transpose() * mu_.head(dim_passive_);
    }
  }
  double error = 0;
  error += q_res_.squaredNorm();
  error += v_res_.squaredNorm();
  error += u_res_.squaredNorm();
  error += lq_.squaredNorm();
  error += lv_.squaredNorm();
  error += la_.squaredNorm();
  error += lu_.squaredNorm();
  error += joint_constraints_.residualSquaredNrom(dtau, q, v, a, u);
  if (dimf > 0) {
    error += lf_.head(dimf).squaredNorm();
  }
  if (dimc > 0) {
    error += C_res_.head(dimc).squaredNorm();
  }
  return error;
}

} // namespace idocp