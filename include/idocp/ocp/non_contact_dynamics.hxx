#ifndef IDOCP_NON_CONTACT_DYNAMICS_HXX_
#define IDOCP_NON_CONTACT_DYNAMICS_HXX_

#include "idocp/ocp/non_contact_dynamics.hpp"

#include <assert.h>

namespace idocp {

inline NonContactDynamics::NonContactDynamics(const Robot& robot) 
  : lu_condensed_(Eigen::VectorXd::Zero(robot.dimv())),
    ID_(Eigen::VectorXd::Zero(robot.dimv())),
    Quu_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dID_dq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dID_dv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dID_da_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Quu_du_dq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Quu_du_dv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Quu_du_da_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dimv_(robot.dimv()) {
}


inline NonContactDynamics::NonContactDynamics() 
  : lu_condensed_(),
    ID_(),
    Quu_(),
    dID_dq_(),
    dID_dv_(),
    dID_da_(),
    Quu_du_dq_(),
    Quu_du_dv_(),
    Quu_du_da_(),
    dimv_(0) {
}


inline NonContactDynamics::~NonContactDynamics() {
}


inline void NonContactDynamics::linearizeNonContactDynamics(
    Robot& robot, const ContactStatus& contact_status, const double dtau, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) { 
  assert(!contact_status.hasActiveContacts());
  assert(dtau > 0);
  linearizeInverseDynamics(robot, contact_status, s, kkt_residual);
  // augment inverse dynamics constraint
  kkt_residual.la.noalias() += dtau * dID_da_.transpose() * s.beta;
  kkt_residual.lq().noalias() += dtau * dID_dq_.transpose() * s.beta;
  kkt_residual.lv().noalias() += dtau * dID_dv_.transpose() * s.beta;
  kkt_residual.lu().noalias() -= dtau * s.beta; 
}


inline void NonContactDynamics::condenseNonContactDynamics(
    Robot& robot, const ContactStatus& contact_status, const double dtau, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  assert(!contact_status.hasActiveContacts());
  assert(dtau > 0);
  lu_condensed_.noalias() = kkt_residual.lu() 
          + kkt_matrix.Quu().diagonal().asDiagonal() * ID_;
  // condense KKT residual
  kkt_residual.la.noalias() = du_da_.transpose() * lu_condensed_;
  kkt_residual.lq().noalias() = du_dq_.transpose() * lu_condensed_;
  kkt_residual.lv().noalias() = du_dv_.transpose() * lu_condensed_;
  kkt_residual.lu().noalias() -= dtau * s.beta;   
  // condense KKT Hessian
  Quu_du_da_.noalias() = kkt_matrix.Quu().diagonal().asDiagonal() * du_da_;
  kkt_matrix.Qaa().noalias() = du_da_.transpose() * Quu_du_da_;
  Quu_du_dq_.noalias() = kkt_matrix.Quu().diagonal().asDiagonal() * du_dq_;
  Quu_du_dv_.noalias() = kkt_matrix.Quu().diagonal().asDiagonal() * du_dv_;
  kkt_matrix.Qaq().noalias() = du_da_.transpose() * Quu_du_dq_;
  kkt_matrix.Qav().noalias() = du_da_.transpose() * Quu_du_dv_;
  kkt_matrix.Qqq().noalias() = du_dq_.transpose() * Quu_du_dq_;
  kkt_matrix.Qqv().noalias() = du_dq_.transpose() * Quu_du_dv_;
  kkt_matrix.Qvv().noalias() = du_dv_.transpose() * Quu_du_dv_;
}


inline void NonContactDynamics::computeCondensedDirection(
    const double dtau, const KKTMatrix& kkt_matrix, 
    const KKTResidual& kkt_residual, SplitDirection& d) {
  assert(dtau > 0);
  d.du() = ID_;
  d.du().noalias() += dID_dq_ * d.dq();
  d.du().noalias() += dID_dv_ * d.dv();
  d.du().noalias() += dID_da_ * d.da();
  d.dbeta().noalias() = (kkt_residual.lu()  + kkt_matrix.Quu() * d.du()) / dtau;
}


inline void NonContactDynamics::computeNonContactDynamicsResidual(
    Robot& robot, const ContactStatus& contact_status, const double dtau, 
    const SplitSolution& s, KKTResidual& kkt_residual) {
  setContactForcesZero(robot, contact_status, s);
  computeInverseDynamicsResidual(robot, s, kkt_residual);
}


inline double NonContactDynamics::l1NormNonContactDynamicsResidual(
    const double dtau) const {
  assert(dtau > 0);
  return (dtau * ID_.lpNorm<1>());
}


inline double NonContactDynamics::squaredNormNonContactDynamicsResidual(
    const double dtau) const {
  assert(dtau > 0);
  return (dtau * dtau * ID_.squaredNorm());
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3, 
          typename MatrixType4>
inline void NonContactDynamics::getStateFeedbackGain(
    const Eigen::MatrixBase<MatrixType1>& Kaq,
    const Eigen::MatrixBase<MatrixType2>& Kav,
    const Eigen::MatrixBase<MatrixType3>& Kuq,
    const Eigen::MatrixBase<MatrixType4>& Kuv) const {
  assert(Kaq.rows() == dimv_);
  assert(Kaq.cols() == dimv_);
  assert(Kav.rows() == dimv_);
  assert(Kav.cols() == dimv_);
  assert(Kuq.rows() == dimv_);
  assert(Kuq.cols() == dimv_);
  assert(Kuv.rows() == dimv_);
  assert(Kuv.cols() == dimv_);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(Kuq) = dID_dq_;
  const_cast<Eigen::MatrixBase<MatrixType3>&>(Kuq).noalias() += dID_da_ * Kaq;
  const_cast<Eigen::MatrixBase<MatrixType4>&>(Kuv) = dID_dv_;
  const_cast<Eigen::MatrixBase<MatrixType4>&>(Kuv).noalias() += dID_da_ * Kav;
}


inline void NonContactDynamics::linearizeInverseDynamics(
    Robot& robot, const ContactStatus& contact_status, const SplitSolution& s, 
    KKTResidual& kkt_residual) {
  setContactForcesZero(robot, contact_status, s);
  computeInverseDynamicsResidual(robot, s, kkt_residual);
  robot.RNEADerivatives(s.q, s.v, s.a, dID_dq_, dID_dv_, dID_da_);
}


inline void NonContactDynamics::setContactForcesZero(
    Robot& robot, const ContactStatus& contact_status, const SplitSolution& s) {
  robot.setContactForces(contact_status, s.f);
}


inline void NonContactDynamics::computeInverseDynamicsResidual(
    Robot& robot, const SplitSolution& s, KKTResidual& kkt_residual) {
  robot.RNEA(s.q, s.v, s.a, ID_);
  ID_.noalias() -= s.u;
}

} // namespace idocp 

#endif // IDOCP_NON_CONTACT_DYNAMICS_HXX_ 