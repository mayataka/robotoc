#ifndef IDOCP_NON_CONTACT_DYNAMICS_HXX_
#define IDOCP_NON_CONTACT_DYNAMICS_HXX_

#include "idocp/ocp/non_contact_dynamics.hpp"

#include <assert.h>

namespace idocp {

inline NonContactDynamics::NonContactDynamics(const Robot& robot) 
  : lu_condensed_(Eigen::VectorXd::Zero(robot.dimv())),
    C_floating_base_(Eigen::VectorXd::Zero(robot.dim_passive())),
    du_dq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    du_dv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    du_da_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Quu_du_dq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Quu_du_dv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Quu_du_da_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    has_floating_base_(robot.has_floating_base()),
    dimv_(robot.dimv()) {
}


inline NonContactDynamics::NonContactDynamics() 
  : lu_condensed_(),
    C_floating_base_(),
    du_dq_(),
    du_dv_(),
    du_da_(),
    Quu_du_dq_(),
    Quu_du_dv_(),
    Quu_du_da_(),
    has_floating_base_(false),
    dimv_(0) {
}


inline NonContactDynamics::~NonContactDynamics() {
}


inline void NonContactDynamics::linearizeRobotDynamics(
    Robot& robot, const ContactStatus& contact_status, const double dtau, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) { 
  assert(!contact_status.hasActiveContacts());
  assert(dtau > 0);
  linearizeInverseDynamics(robot, contact_status, s, kkt_residual);
  // augment inverse dynamics constraint
  kkt_residual.la().noalias() += dtau * dID_da_.transpose() * s.beta;
  kkt_residual.lq().noalias() += dtau * dID_dq_.transpose() * s.beta;
  kkt_residual.lv().noalias() += dtau * dID_dv_.transpose() * s.beta;
  kkt_residual.lu.noalias() -= dtau * s.beta; 
  // // augment floating base constraint
  // if (has_floating_base_) {
  //   computeFloatingBaseConstraintResidual(robot, dtau, s, kkt_residual);
  //   kkt_residual.lu.template head<kDimFloatingBase>().noalias() 
  //       += dtau * s.mu_floating_base();
  // }
}


inline void NonContactDynamics::condenseRobotDynamics(
    Robot& robot, const ContactStatus& contact_status, const double dtau, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  assert(!contact_status.hasActiveContacts());
  assert(dtau > 0);
  linearizeInverseDynamics(robot, contact_status, s, kkt_residual);
  lu_condensed_.noalias() 
      = kkt_residual.lu 
          + kkt_matrix.Quu.diagonal().asDiagonal() * kkt_residual.u_res;
  // condense Newton residual
  kkt_residual.la().noalias() = du_da_.transpose() * lu_condensed_;
  kkt_residual.lq().noalias() = du_dq_.transpose() * lu_condensed_;
  kkt_residual.lv().noalias() = du_dv_.transpose() * lu_condensed_;
  kkt_residual.lu.noalias() -= dtau * s.beta;   
  // condense Hessian
  // Assume that Quu has only diagonal elements, hence we write as 
  // kkt_matrix.Quu.dianobal().asDiagonal() for efficiency.
  // Otherwise, simply use kkt_matrix.Quu instead.
  Quu_du_da_.noalias() = kkt_matrix.Quu.diagonal().asDiagonal() * du_da_;
  kkt_matrix.Qaa().noalias() = du_da_.transpose() * Quu_du_da_;
  Quu_du_dq_.noalias() = kkt_matrix.Quu.diagonal().asDiagonal() * du_dq_;
  Quu_du_dv_.noalias() = kkt_matrix.Quu.diagonal().asDiagonal() * du_dv_;
  kkt_matrix.Qaq().noalias() = du_da_.transpose() * Quu_du_dq_;
  kkt_matrix.Qav().noalias() = du_da_.transpose() * Quu_du_dv_;
  kkt_matrix.Qqq().noalias() = du_dq_.transpose() * Quu_du_dq_;
  kkt_matrix.Qqv().noalias() = du_dq_.transpose() * Quu_du_dv_;
  kkt_matrix.Qvv().noalias() = du_dv_.transpose() * Quu_du_dv_;
  // // condense floating base constraint
  // if (has_floating_base_) {
  //   kkt_residual.C_floating_base().noalias() 
  //       = dtau * (kkt_residual.u_res.template head<kDimFloatingBase>()
  //                 + s.u.template head<kDimFloatingBase>());
  //   // Not condensed constraint residual to use in l1NormInverseDynamicsResidual
  //   // and squaredNormNonContactDynamicsResidual.
  //   C_floating_base_ = dtau * s.u.template head<kDimFloatingBase>();
  //   kkt_matrix.Ca_floating_base()
  //       = dtau * du_da_.template topRows<kDimFloatingBase>();
  //   kkt_matrix.Cq_floating_base() 
  //       = dtau * du_dq_.template topRows<kDimFloatingBase>();
  //   kkt_matrix.Cv_floating_base() 
  //       = dtau * du_dv_.template topRows<kDimFloatingBase>();
  //   kkt_residual.la().noalias() += kkt_matrix.Ca().transpose() * s.mu_stack();
  //   kkt_residual.lq().noalias() += kkt_matrix.Cq().transpose() * s.mu_stack();
  //   kkt_residual.lv().noalias() += kkt_matrix.Cv().transpose() * s.mu_stack();
  // }
}


inline void NonContactDynamics::computeCondensedDirection(
    const double dtau, const KKTMatrix& kkt_matrix, 
    const KKTResidual& kkt_residual, SplitDirection& d) {
  assert(dtau > 0);
  d.du = kkt_residual.u_res;
  d.du.noalias() += du_dq_ * d.dq();
  d.du.noalias() += du_dv_ * d.dv();
  d.du.noalias() += du_da_ * d.da();
  d.dbeta.noalias() = (kkt_residual.lu  + kkt_matrix.Quu * d.du) / dtau;
  // if (has_floating_base_) {
  //   d.dbeta.template head<kDimFloatingBase>().noalias() 
  //       += d.dmu().template head<kDimFloatingBase>();
  // }
}


inline void NonContactDynamics::computeRobotDynamicsResidual(
    Robot& robot, const ContactStatus& contact_status, const double dtau, 
    const SplitSolution& s, KKTResidual& kkt_residual) {
  setContactForcesZero(robot, contact_status, s);
  computeInverseDynamicsResidual(robot, s, kkt_residual);
  computeFloatingBaseConstraintResidual(robot, dtau, s, kkt_residual);
}


inline void NonContactDynamics::linearizeInverseDynamics(
    Robot& robot, const ContactStatus& contact_status, const SplitSolution& s, 
    KKTResidual& kkt_residual) {
  setContactForcesZero(robot, contact_status, s);
  computeInverseDynamicsResidual(robot, s, kkt_residual);
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq_, du_dv_, du_da_);
}


inline void NonContactDynamics::setContactForcesZero(
    Robot& robot, const ContactStatus& contact_status, const SplitSolution& s) {
  robot.setContactForces(contact_status, s.f);
}


inline void NonContactDynamics::computeInverseDynamicsResidual(
    Robot& robot, const SplitSolution& s, KKTResidual& kkt_residual) {
  robot.RNEA(s.q, s.v, s.a, kkt_residual.u_res);
  kkt_residual.u_res.noalias() -= s.u;
}


inline void NonContactDynamics::computeFloatingBaseConstraintResidual(
    const Robot& robot, const double dtau, 
    const SplitSolution& s, KKTResidual& kkt_residual) {
  if (has_floating_base_) {
    C_floating_base_ = dtau * s.u.template head<kDimFloatingBase>();
    kkt_residual.C_floating_base() = C_floating_base_;
  }
}


inline double NonContactDynamics::l1NormNonContactDynamicsResidual(
    const double dtau, const KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  return (dtau * kkt_residual.u_res.lpNorm<1>() 
          + C_floating_base_.lpNorm<1>());
}


inline double NonContactDynamics::squaredNormNonContactDynamicsResidual(
    const double dtau, const KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  return (dtau * dtau * kkt_residual.u_res.squaredNorm() 
          + C_floating_base_.squaredNorm());
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3, 
          typename MatrixType4>
inline void NonContactDynamics::getStateFeedbackGain(
    const Eigen::MatrixBase<MatrixType1>& da_dq,
    const Eigen::MatrixBase<MatrixType2>& da_dv,
    const Eigen::MatrixBase<MatrixType3>& Kuq,
    const Eigen::MatrixBase<MatrixType4>& Kuv) const {
  assert(da_dq.rows() == da_dq.cols());
  assert(da_dv.rows() == da_dv.cols());
  assert(df_dq.rows() == dimf_);
  assert(df_dv.rows() == dimf_);
  assert(Kuq.rows() == Kuq.cols());
  assert(Kuv.rows() == Kuv.cols());
  assert(da_dq.rows() == da_dv.rows());
  assert(da_dq.rows() == df_dq.cols());
  assert(da_dq.rows() == df_dv.cols());
  assert(da_dq.rows() == Kuq.rows());
  assert(da_dq.rows() == Kuv.rows());
  const_cast<Eigen::MatrixBase<MatrixType3>&>(Kuq) = du_dq_;
  const_cast<Eigen::MatrixBase<MatrixType3>&>(Kuq).noalias() += du_da_ * da_dq;
  const_cast<Eigen::MatrixBase<MatrixType4>&>(Kuv) = du_dv_;
  const_cast<Eigen::MatrixBase<MatrixType4>&>(Kuv).noalias() += du_da_ * da_dv;
}

} // namespace idocp 

#endif // IDOCP_NON_CONTACT_DYNAMICS_HXX_ 