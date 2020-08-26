#ifndef IDOCP_ROBOT_DYNAMICS_HXX_
#define IDOCP_ROBOT_DYNAMICS_HXX_

#include <assert.h>

namespace idocp {

inline RobotDynamics::RobotDynamics(const Robot& robot) 
  : lu_condensed_(Eigen::VectorXd::Zero(robot.dimv())),
    du_dq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    du_dv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    du_da_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    du_df_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    has_floating_base_(robot.has_floating_base()),
    has_active_contacts_(robot.has_active_contacts()),
    dimf_(robot.dimf()) {
}


inline RobotDynamics::RobotDynamics() 
  : lu_condensed_(),
    du_dq_(),
    du_dv_(),
    du_da_(),
    du_df_(),
    has_floating_base_(false),
    has_active_contacts_(false),
    dimf_(0) {
}


inline RobotDynamics::~RobotDynamics() {
}


inline void RobotDynamics::augmentRobotDynamics(Robot& robot, const double dtau, 
                                                const SplitSolution& s, 
                                                KKTMatrix& kkt_matrix, 
                                                KKTResidual& kkt_residual) {
  setContactStatus(robot);
  linearizeInverseDynamics(robot, s, kkt_residual);
  // augment inverse dynamics
  kkt_residual.la().noalias() += dtau * du_da_.transpose() * s.beta;
  if (robot.has_active_contacts()) {
    kkt_residual.lf().noalias() += dtau * du_df_active_().transpose() * s.beta;
  }
  kkt_residual.lq().noalias() += dtau * du_dq_.transpose() * s.beta;
  kkt_residual.lv().noalias() += dtau * du_dv_.transpose() * s.beta;
  kkt_residual.lu.noalias() -= dtau * s.beta; 
  // augment floating base constraint
  linearizeFloatingBaseConstraint(robot, dtau, s, kkt_residual);
  // augment contact constraint
  if (robot.has_active_contacts()) {
    linearizeContactConstraint(robot, dtau, kkt_matrix, kkt_residual);
    kkt_residual.la().noalias() 
        += kkt_matrix.Ca().topRows(robot.dimf()).transpose() 
            * s.mu.head(robot.dimf());
    kkt_residual.lq().noalias()
        += kkt_matrix.Cq().topRows(robot.dimf()).transpose() 
            * s.mu.head(robot.dimf());
    kkt_residual.lv().noalias() 
        += kkt_matrix.Cv().topRows(robot.dimf()).transpose() 
            * s.mu.head(robot.dimf());
  }
}


inline void RobotDynamics::condenseRobotDynamics(Robot& robot, 
                                                 const double dtau, 
                                                 const SplitSolution& s, 
                                                 KKTMatrix& kkt_matrix, 
                                                 KKTResidual& kkt_residual) {
  setContactStatus(robot);
  linearizeInverseDynamics(robot, s, kkt_residual);
  if (robot.has_floating_base()) {
    kkt_residual.lu.template head<kDimFloatingBase>().noalias()
        += dtau * s.mu_active().template tail<kDimFloatingBase>();
  }
  lu_condensed_ = kkt_residual.lu + kkt_matrix.Quu * kkt_residual.u_res;
  // condense Newton residual
  kkt_residual.la() = du_da_.transpose() * lu_condensed_;
  if (robot.has_active_contacts()) {
    kkt_residual.lf() = du_df_active_().transpose() * lu_condensed_;
  }
  kkt_residual.lq() = du_dq_.transpose() * lu_condensed_;
  kkt_residual.lv() = du_dv_.transpose() * lu_condensed_;
  kkt_residual.lu.noalias() -= dtau * s.beta;   
  // condense Hessian
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
  // condense floating base constraint
  if (robot.has_floating_base()) {
    kkt_residual.C().template tail<kDimFloatingBase>() 
        = dtau * (kkt_residual.u_res.template head<kDimFloatingBase>()
                  + s.u.template head<kDimFloatingBase>());
    kkt_matrix.Ca().template bottomRows<kDimFloatingBase>() 
        = dtau * du_da_.template topRows<kDimFloatingBase>();
    if (robot.has_active_contacts()) {
      kkt_matrix.Cf().template bottomRows<kDimFloatingBase>() 
          = dtau * du_df_active_().template topRows<kDimFloatingBase>();
    }
    kkt_matrix.Cq().template bottomRows<kDimFloatingBase>() 
        = dtau * du_dq_.template topRows<kDimFloatingBase>();
    kkt_matrix.Cv().template bottomRows<kDimFloatingBase>() 
        = dtau * du_dv_.template topRows<kDimFloatingBase>();
  }
  if (robot.has_active_contacts()) {
    linearizeContactConstraint(robot, dtau, kkt_matrix, kkt_residual);
  }
  // augment floating base constraint and contact constraint together
  if (robot.has_floating_base() || robot.has_active_contacts()) {
    kkt_residual.la().noalias() += kkt_matrix.Ca().transpose() * s.mu_active();
    if (robot.has_floating_base()) {
      kkt_residual.lf().noalias() 
          += kkt_matrix.Cf().template bottomRows<kDimFloatingBase>().transpose() 
              * s.mu_active().template tail<kDimFloatingBase>();
    }
    kkt_residual.lq().noalias() += kkt_matrix.Cq().transpose() * s.mu_active();
    kkt_residual.lv().noalias() += kkt_matrix.Cv().transpose() * s.mu_active();
  }
}


inline void RobotDynamics::computeCondensedDirection(
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
    d.dbeta.template head<kDimFloatingBase>().noalias() 
        += d.dmu().template tail<kDimFloatingBase>();
  }
}


inline double RobotDynamics::violationL1Norm(
    Robot& robot, const double dtau, const SplitSolution& s, 
    KKTResidual& kkt_residual) const {
  double violation = dtau * kkt_residual.u_res.lpNorm<1>();
  if (has_floating_base_) {
    violation += dtau * s.u.template head<kDimFloatingBase>().lpNorm<1>();
  }
  if (has_active_contacts_) {
    violation += kkt_residual.C().head(dimf_).lpNorm<1>();
  }
  return violation;
}


inline double RobotDynamics::computeViolationL1Norm(
    Robot& robot, const double dtau, const SplitSolution& s, 
    KKTResidual& kkt_residual) const {
  if (robot.has_active_contacts()) {
    robot.setContactForces(s.f);
  }
  robot.RNEA(s.q, s.v, s.a, kkt_residual.u_res);
  kkt_residual.u_res.noalias() -= s.u;
  double violation = dtau * kkt_residual.u_res.lpNorm<1>();
  if (robot.has_floating_base()) {
    violation += dtau * s.u.template head<kDimFloatingBase>().lpNorm<1>();
  }
  if (robot.has_active_contacts()) {
    robot.computeBaumgarteResidual(dtau, kkt_residual.C());
    violation += kkt_residual.C().head(robot.dimf()).lpNorm<1>();
  }
  return violation;
}


inline void RobotDynamics::setContactStatus(const Robot& robot) {
  dimf_ = robot.dimf();
  has_active_contacts_ = robot.has_active_contacts();
}


inline void RobotDynamics::linearizeInverseDynamics(Robot& robot, 
                                                    const SplitSolution& s, 
                                                    KKTResidual& kkt_residual) {
  if (robot.has_active_contacts()) {
    robot.setContactForces(s.f);
  }
  robot.RNEA(s.q, s.v, s.a, kkt_residual.u_res);
  kkt_residual.u_res.noalias() -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq_, du_dv_, du_da_);
  if (robot.has_active_contacts()) {
    robot.dRNEAPartialdFext(du_df_);
  }
}


inline void RobotDynamics::linearizeContactConstraint(
    Robot& robot, const double dtau, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  robot.computeBaumgarteResidual(dtau, kkt_residual.C());
  robot.computeBaumgarteDerivatives(dtau, kkt_matrix.Cq(), kkt_matrix.Cv(), 
                                    kkt_matrix.Ca());
}


inline void RobotDynamics::linearizeFloatingBaseConstraint(
    const Robot& robot, const double dtau, const SplitSolution& s, 
    KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  if (robot.has_floating_base()) {
    kkt_residual.C().template tail<kDimFloatingBase>() 
        = dtau * s.u.template head<kDimFloatingBase>();
    kkt_residual.lu.template head<kDimFloatingBase>().noalias()
          += dtau * s.mu_active().template tail<kDimFloatingBase>();
  }
}


inline Eigen::Block<Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic, true> 
RobotDynamics::du_df_active_() {
  return du_df_.leftCols(dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd, Eigen::Dynamic, 
                          Eigen::Dynamic, true> 
RobotDynamics::du_df_active_() const {
  return du_df_.leftCols(dimf_);
}

} // namespace idocp 

#endif // IDOCP_ROBOT_DYNAMICS_HXX_