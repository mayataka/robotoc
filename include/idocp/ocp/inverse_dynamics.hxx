#ifndef IDOCP_INVERSE_DYNAMICS_HXX_
#define IDOCP_INVERSE_DYNAMICS_HXX_

namespace idocp {

inline InverseDynamics::InverseDynamics(const Robot& robot) 
  : lu_condensed_(Eigen::VectorXd::Zero(robot.dimv())),
    du_dq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    du_dv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    du_da_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    du_df_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    has_floating_base_(robot.has_floating_base()),
    dimv_(robot.dimv()),
    dim_passive_(robot.dim_passive()),
    dimf_(robot.dimf()) {
}


inline InverseDynamics::InverseDynamics() 
  : lu_condensed_(),
    du_dq_(),
    du_dv_(),
    du_da_(),
    du_df_(),
    has_floating_base_(false),
    dimv_(0),
    dim_passive_(0),
    dimf_(0) {
}


inline InverseDynamics::~InverseDynamics() {
}


inline void InverseDynamics::augmentInverseDynamics(Robot& robot, 
                                                    const double dtau, 
                                                    const SplitSolution& s, 
                                                    KKTResidual& kkt_residual) {
  setContactStatus(robot);
  if (robot.has_active_contacts()) {
    robot.setContactForces(s.f);
  }
  robot.RNEA(s.q, s.v, s.a, kkt_residual.u_res);
  kkt_residual.u_res.noalias() -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq_, du_dv_, du_da_);
  if (robot.has_active_contacts()) {
    robot.dRNEAPartialdFext(du_df_);
  }
  kkt_residual.lq().noalias() += dtau * du_dq_.transpose() * s.beta;
  kkt_residual.lv().noalias() += dtau * du_dv_.transpose() * s.beta;
  kkt_residual.la().noalias() += dtau * du_da_.transpose() * s.beta;
  if (robot.has_active_contacts()) {
    kkt_residual.lf().noalias() += dtau * du_df_active_().transpose() * s.beta;
  }
  kkt_residual.lu.noalias() -= dtau * s.beta; 
}


inline void InverseDynamics::condenseInverseDynamics(
    Robot& robot, const double dtau, const SplitSolution& s, 
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  setContactStatus(robot);
  if (robot.has_active_contacts()) {
    robot.setContactForces(s.f);
  }
  robot.RNEA(s.q, s.v, s.a, kkt_residual.u_res);
  kkt_residual.u_res.noalias() -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq_, du_dv_, du_da_);
  if (robot.has_active_contacts()) {
    robot.dRNEAPartialdFext(du_df_);
  }
  if (robot.has_floating_base()) {
    kkt_residual.C().tail(robot.dim_passive()) 
        = dtau * s.u.head(robot.dim_passive());
    kkt_residual.lu.head(robot.dim_passive()).noalias()
        += dtau * s.mu_active().tail(robot.dim_passive());
  }
  lu_condensed_ = kkt_residual.lu + kkt_matrix.Quu * kkt_residual.u_res;
  // condense Newton residual
  kkt_residual.lu.noalias() -= dtau * s.beta;   
  kkt_residual.lq() = du_dq_.transpose() * lu_condensed_;
  kkt_residual.lv() = du_dv_.transpose() * lu_condensed_;
  kkt_residual.la() = du_da_.transpose() * lu_condensed_;
  if (robot.has_active_contacts()) {
    kkt_residual.lf() = du_df_active_().transpose() * lu_condensed_;
  }
  kkt_matrix.Qaa() = du_da_.transpose() * kkt_matrix.Quu * du_da_;
  if (robot.has_active_contacts()) {
    kkt_matrix.Qaf() = du_da_.transpose() * kkt_matrix.Quu * du_df_active_(); 
  }
  kkt_matrix.Qaq() = du_da_.transpose() * kkt_matrix.Quu * du_dq_;
  kkt_matrix.Qav() = du_da_.transpose() * kkt_matrix.Quu * du_dv_;
  if (robot.has_active_contacts()) {
    kkt_matrix.Qff() 
        = du_df_active_().transpose() * kkt_matrix.Quu * du_df_active_();
    kkt_matrix.Qfq() = du_df_active_().transpose() * kkt_matrix.Quu * du_dq_;
    kkt_matrix.Qfv() = du_df_active_().transpose() * kkt_matrix.Quu * du_dv_;
  }
  kkt_matrix.Qqq() = du_dq_.transpose() * kkt_matrix.Quu * du_dq_;
  kkt_matrix.Qqv() = du_dq_.transpose() * kkt_matrix.Quu * du_dv_;
  kkt_matrix.Qvv() = du_dv_.transpose() * kkt_matrix.Quu * du_dv_;
  if (robot.has_floating_base()) {
    kkt_residual.C().tail(dim_passive_).noalias()
        += dtau * kkt_residual.u_res.head(dim_passive_);
    kkt_matrix.Cq().bottomRows(dim_passive_) 
        = dtau * du_dq_.topRows(dim_passive_);
    kkt_matrix.Cv().bottomRows(dim_passive_) 
        = dtau * du_dv_.topRows(dim_passive_);
    kkt_matrix.Ca().bottomRows(dim_passive_) 
        = dtau * du_da_.topRows(dim_passive_);
    kkt_matrix.Cf().bottomRows(dim_passive_) 
        = dtau * du_df_active_().topRows(dim_passive_);
  }
}


inline void InverseDynamics::computeCondensedDirection(
    const double dtau, const KKTMatrix& kkt_matrix, 
    const KKTResidual& kkt_residual, SplitDirection& d) {
  assert(dtau > 0);
  d.du = kkt_residual.u_res;
  d.du.noalias() += du_dq_ * d.dq();
  d.du.noalias() += du_dv_ * d.dv();
  d.du.noalias() += du_da_ * d.da();
  if (has_active_contacts_) {
    d.du.noalias() += du_df_active_() * d.df();
  }
  d.dbeta = (kkt_residual.lu  + kkt_matrix.Quu * d.du) / dtau;
  if (has_floating_base_) {
    d.dbeta.head(dim_passive_).noalias() += d.dmu().tail(dim_passive_);
  }
}


inline double InverseDynamics::violationL1Norm(
    const double dtau, const KKTResidual& kkt_residual) const {
  return dtau * kkt_residual.u_res.lpNorm<1>();
}


inline double InverseDynamics::computeViolationL1Norm(
    Robot& robot, const double dtau, const SplitSolution& s, 
    KKTResidual& kkt_residual) const {
  if (robot.has_active_contacts()) {
    robot.setContactForces(s.f);
  }
  robot.RNEA(s.q, s.v, s.a, kkt_residual.u_res);
  kkt_residual.u_res.noalias() -= s.u;
  return dtau * kkt_residual.u_res.lpNorm<1>();
}


inline void InverseDynamics::setContactStatus(const Robot& robot) {
  has_active_contacts_ = robot.has_active_contacts();
  dimf_ = robot.dimf();
}


inline Eigen::Block<Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic, true> 
InverseDynamics::du_df_active_() {
  return du_df_.leftCols(dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd, Eigen::Dynamic, 
                          Eigen::Dynamic, true> 
InverseDynamics::du_df_active_() const {
  return du_df_.leftCols(dimf_);
}

} // namespace idocp

#endif // IDOCP_INVERSE_DYNAMICS_HXX_ 