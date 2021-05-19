#ifndef IDOCP_UNCONSTRAINED_DYNAMICS_HXX_
#define IDOCP_UNCONSTRAINED_DYNAMICS_HXX_

#include "idocp/unocp/unconstrained_dynamics.hpp"

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
    if (robot.hasFloatingBase()) {
      throw std::logic_error("robot has floating base: robot should have no constraints!");
    }
    if (robot.maxPointContacts() > 0) {
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
    Robot& robot, const double dt, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) { 
  assert(dt > 0);
  linearizeInverseDynamics(robot, s);
  // augment inverse dynamics constraint
  kkt_residual.lq().noalias() += dt * dID_dq_.transpose() * s.beta;
  kkt_residual.lv().noalias() += dt * dID_dv_.transpose() * s.beta;
  kkt_residual.la.noalias()   += dt * dID_da_.transpose() * s.beta;
  kkt_residual.lu.noalias() -= dt * s.beta; 
}


inline void UnconstrainedDynamics::condenseUnconstrainedDynamics(
    const SplitKKTMatrix& kkt_matrix, const SplitKKTResidual& kkt_residual, 
    SplitUnKKTMatrix& unkkt_matrix, SplitUnKKTResidual& unkkt_residual) {
  // condense KKT residual
  lu_condensed_ = kkt_residual.lu;
  lu_condensed_.noalias() += kkt_matrix.Quu.diagonal().asDiagonal() * ID_;
  unkkt_residual.lq() = kkt_residual.lq();
  unkkt_residual.lv() = kkt_residual.lv();
  unkkt_residual.la() = kkt_residual.la; 
  unkkt_residual.lq().noalias() += dID_dq_.transpose() * lu_condensed_;
  unkkt_residual.lv().noalias() += dID_dv_.transpose() * lu_condensed_;
  unkkt_residual.la().noalias() += dID_da_.transpose() * lu_condensed_;
  unkkt_residual.Fx() = kkt_residual.Fx;
  // condense KKT Hessian
  Quu_dID_dq_.noalias() = kkt_matrix.Quu.diagonal().asDiagonal() * dID_dq_;
  Quu_dID_dv_.noalias() = kkt_matrix.Quu.diagonal().asDiagonal() * dID_dv_;
  Quu_dID_da_.noalias() = kkt_matrix.Quu.diagonal().asDiagonal() * dID_da_;
  unkkt_matrix.Qqq().noalias() = dID_dq_.transpose() * Quu_dID_dq_;
  unkkt_matrix.Qqv().noalias() = dID_dq_.transpose() * Quu_dID_dv_;
  unkkt_matrix.Qvv().noalias() = dID_dv_.transpose() * Quu_dID_dv_;
  unkkt_matrix.Qaq().noalias() = dID_da_.transpose() * Quu_dID_dq_;
  unkkt_matrix.Qav().noalias() = dID_da_.transpose() * Quu_dID_dv_;
  unkkt_matrix.Qaa().noalias() = dID_da_.transpose() * Quu_dID_da_;
  unkkt_matrix.Qqq().noalias() += kkt_matrix.Qqq();
  unkkt_matrix.Qvv().diagonal().noalias() += kkt_matrix.Qvv().diagonal();
  unkkt_matrix.Qaa().diagonal().noalias() += kkt_matrix.Qaa.diagonal();
}


inline void UnconstrainedDynamics::computeCondensedDirection(
    const double dt, const SplitKKTMatrix& kkt_matrix, 
    const SplitKKTResidual& kkt_residual, SplitDirection& d) {
  assert(dt > 0);
  d.du = ID_;
  d.du.noalias() += dID_dq_ * d.dq();
  d.du.noalias() += dID_dv_ * d.dv();
  d.du.noalias() += dID_da_ * d.da();
  d.dbeta().noalias() = (kkt_residual.lu  + kkt_matrix.Quu * d.du) / dt;
}


inline void UnconstrainedDynamics::computeUnconstrainedDynamicsResidual(
    Robot& robot, const SplitSolution& s) {
  computeInverseDynamicsResidual(robot, s);
}


inline double UnconstrainedDynamics::l1NormUnconstrainedDynamicsResidual(
    const double dt) const {
  assert(dt > 0);
  return (dt * ID_.lpNorm<1>());
}


inline double UnconstrainedDynamics::squaredNormUnconstrainedDynamicsResidual(
    const double dt) const {
  assert(dt > 0);
  return (dt * dt * ID_.squaredNorm());
}


template <typename MatrixType1, typename MatrixType2>
inline void UnconstrainedDynamics::getStateFeedbackGain(
    const Eigen::MatrixBase<MatrixType1>& Ka,
    const Eigen::MatrixBase<MatrixType2>& Ku) const {
  assert(Ka.rows() == dimv_);
  assert(Ka.cols() == 2*dimv_);
  assert(Ku.rows() == dimv_);
  assert(Ku.cols() == 2*dimv_);
  getStateFeedbackGain(
      Ka.leftCols(dimv_), Ka.rightCols(dimv_),
      const_cast<Eigen::MatrixBase<MatrixType2>&>(Ku).leftCols(dimv_),
      const_cast<Eigen::MatrixBase<MatrixType2>&>(Ku).rightCols(dimv_));
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