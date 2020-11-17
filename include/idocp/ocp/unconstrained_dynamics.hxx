#ifndef IDOCP_UNCONSTRAINED_DYNAMICS_HXX_
#define IDOCP_UNCONSTRAINED_DYNAMICS_HXX_

#include "idocp/ocp/unconstrained_dynamics.hpp"

#include <stdexcept>
#include <cassert>

namespace idocp {

inline UnconstrainedDynamics::UnconstrainedDynamics(const Robot& robot) 
  : lu_condensed_(Eigen::VectorXd::Zero(robot.dimv())),
    ID_(Eigen::VectorXd::Zero(robot.dimv())),
    Quu_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dID_dq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dID_dv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dID_da_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Quu_dID_dq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Quu_dID_dv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Quu_dID_da_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dimv_(robot.dimv()) {
  try {
    if (robot.has_floating_base()) {
      throw std::logic_error("robot has floating base: robot should have no constraints!");
    }
    if (robot.max_point_contacts() > 0) {
      throw std::logic_error("robot can have contacts: robot should have no constraints!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline UnconstrainedDynamics::UnconstrainedDynamics() 
  : lu_condensed_(),
    ID_(),
    Quu_(),
    dID_dq_(),
    dID_dv_(),
    dID_da_(),
    Quu_dID_dq_(),
    Quu_dID_dv_(),
    Quu_dID_da_(),
    dimv_(0) {
}


inline UnconstrainedDynamics::~UnconstrainedDynamics() {
}


inline void UnconstrainedDynamics::linearizeUnconstrainedDynamics(
    Robot& robot, const double dtau, const SplitSolution& s, 
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) { 
  assert(dtau > 0);
  linearizeInverseDynamics(robot, s);
  // augment inverse dynamics constraint
  kkt_residual.la.noalias() += dtau * dID_da_.transpose() * s.beta;
  kkt_residual.lq().noalias() += dtau * dID_dq_.transpose() * s.beta;
  kkt_residual.lv().noalias() += dtau * dID_dv_.transpose() * s.beta;
  kkt_residual.lu().noalias() -= dtau * s.beta; 
}


inline void UnconstrainedDynamics::condenseUnconstrainedDynamics(
    Robot& robot, const double dtau, const SplitSolution& s, 
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  assert(dtau > 0);
  lu_condensed_.noalias() = kkt_residual.lu() 
          + kkt_matrix.Quu().diagonal().asDiagonal() * ID_;
  // condense KKT residual
  kkt_residual.la.noalias() += dID_da_.transpose() * lu_condensed_;
  kkt_residual.lq().noalias() += dID_dq_.transpose() * lu_condensed_;
  kkt_residual.lv().noalias() += dID_dv_.transpose() * lu_condensed_;
  kkt_residual.lu().noalias() -= dtau * s.beta;   
  // condense KKT Hessian
  Quu_dID_da_.noalias() = kkt_matrix.Quu().diagonal().asDiagonal() * dID_da_;
  kkt_matrix.Qaa().noalias() += dID_da_.transpose() * Quu_dID_da_;
  Quu_dID_dq_.noalias() = kkt_matrix.Quu().diagonal().asDiagonal() * dID_dq_;
  Quu_dID_dv_.noalias() = kkt_matrix.Quu().diagonal().asDiagonal() * dID_dv_;
  kkt_matrix.Qaq().noalias() += dID_da_.transpose() * Quu_dID_dq_;
  kkt_matrix.Qav().noalias() += dID_da_.transpose() * Quu_dID_dv_;
  kkt_matrix.Qqq().noalias() += dID_dq_.transpose() * Quu_dID_dq_;
  kkt_matrix.Qqv().noalias() += dID_dq_.transpose() * Quu_dID_dv_;
  kkt_matrix.Qvv().noalias() += dID_dv_.transpose() * Quu_dID_dv_;
}


inline void UnconstrainedDynamics::computeCondensedDirection(
    const double dtau, const KKTMatrix& kkt_matrix, 
    const KKTResidual& kkt_residual, SplitDirection& d) {
  assert(dtau > 0);
  d.du() = ID_;
  d.du().noalias() += dID_dq_ * d.dq();
  d.du().noalias() += dID_dv_ * d.dv();
  d.du().noalias() += dID_da_ * d.da();
  d.dbeta().noalias() = (kkt_residual.lu()  + kkt_matrix.Quu() * d.du()) / dtau;
}


inline void UnconstrainedDynamics::computeUnconstrainedDynamicsResidual(
    Robot& robot, const SplitSolution& s) {
  computeInverseDynamicsResidual(robot, s);
}


inline double UnconstrainedDynamics::l1NormUnconstrainedDynamicsResidual(
    const double dtau) const {
  assert(dtau > 0);
  return (dtau * ID_.lpNorm<1>());
}


inline double UnconstrainedDynamics::squaredNormUnconstrainedDynamicsResidual(
    const double dtau) const {
  assert(dtau > 0);
  return (dtau * dtau * ID_.squaredNorm());
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3, 
          typename MatrixType4>
inline void UnconstrainedDynamics::getStateFeedbackGain(
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


inline void UnconstrainedDynamics::linearizeInverseDynamics(
    Robot& robot, const SplitSolution& s) {
  computeInverseDynamicsResidual(robot, s);
  robot.RNEADerivatives(s.q, s.v, s.a, dID_dq_, dID_dv_, dID_da_);
}


inline void UnconstrainedDynamics::computeInverseDynamicsResidual(
    Robot& robot, const SplitSolution& s) {
  robot.RNEA(s.q, s.v, s.a, ID_);
  ID_.noalias() -= s.u;
}

} // namespace idocp 

#endif // IDOCP_UNCONSTRAINED_DYNAMICS_HXX_ 