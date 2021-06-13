#ifndef IDOCP_UNCONSTR_DYNAMICS_HXX_
#define IDOCP_UNCONSTR_DYNAMICS_HXX_

#include "idocp/unconstr/unconstr_dynamics.hpp"

#include <stdexcept>
#include <cassert>

namespace idocp {

inline UnconstrDynamics::UnconstrDynamics(const Robot& robot) 
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


inline UnconstrDynamics::UnconstrDynamics() 
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


inline UnconstrDynamics::~UnconstrDynamics() {
}


inline void UnconstrDynamics::computeUnconstrDynamicsResidual(
    Robot& robot, const SplitSolution& s) {
  robot.RNEA(s.q, s.v, s.a, ID_);
  ID_.noalias() -= s.u;
}


inline void UnconstrDynamics::linearizeUnconstrDynamics(
    Robot& robot, const double dt, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) { 
  assert(dt > 0);
  computeUnconstrDynamicsResidual(robot, s);
  // augment inverse dynamics constraint
  robot.RNEADerivatives(s.q, s.v, s.a, dID_dq_, dID_dv_, dID_da_);
  kkt_residual.lq().noalias() += dt * dID_dq_.transpose() * s.beta;
  kkt_residual.lv().noalias() += dt * dID_dv_.transpose() * s.beta;
  kkt_residual.la.noalias()   += dt * dID_da_.transpose() * s.beta;
  kkt_residual.lu.noalias()   -= dt * s.beta; 
}


inline void UnconstrDynamics::condenseUnconstrDynamics(
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) {
  // condense KKT residual
  lu_condensed_ = kkt_residual.lu;
  lu_condensed_.noalias() += kkt_matrix.Quu.diagonal().asDiagonal() * ID_;
  kkt_residual.lq().noalias() += dID_dq_.transpose() * lu_condensed_;
  kkt_residual.lv().noalias() += dID_dv_.transpose() * lu_condensed_;
  kkt_residual.la.noalias()   += dID_da_.transpose() * lu_condensed_;
  // condense KKT Hessian
  Quu_dID_dq_.noalias() = kkt_matrix.Quu.diagonal().asDiagonal() * dID_dq_;
  Quu_dID_dv_.noalias() = kkt_matrix.Quu.diagonal().asDiagonal() * dID_dv_;
  Quu_dID_da_.noalias() = kkt_matrix.Quu.diagonal().asDiagonal() * dID_da_;
  kkt_matrix.Qqq().noalias() += dID_dq_.transpose() * Quu_dID_dq_;
  kkt_matrix.Qqv().noalias() += dID_dq_.transpose() * Quu_dID_dv_;
  kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
  kkt_matrix.Qvv().noalias() += dID_dv_.transpose() * Quu_dID_dv_;
  kkt_matrix.Qaa.noalias()   += dID_da_.transpose() * Quu_dID_da_;
  //  this is actually Qqa and Qqv
  kkt_matrix.Qqu().noalias() = Quu_dID_dq_.transpose() * dID_da_;
  kkt_matrix.Qvu().noalias() = Quu_dID_dv_.transpose() * dID_da_;
}


inline void UnconstrDynamics::expandPrimal(SplitDirection& d) const {
  d.du = ID_;
  d.du.noalias() += dID_dq_ * d.dq();
  d.du.noalias() += dID_dv_ * d.dv();
  d.du.noalias() += dID_da_ * d.da();
}



inline void UnconstrDynamics::expandDual(const double dt, 
                                         const SplitKKTMatrix& kkt_matrix, 
                                         const SplitKKTResidual& kkt_residual, 
                                         SplitDirection& d) {
  assert(dt > 0);
  d.dbeta().noalias() = (kkt_residual.lu  + kkt_matrix.Quu * d.du) / dt;
}


inline double UnconstrDynamics::squaredNormKKTResidual() const {
  return ID_.squaredNorm();
}


inline double UnconstrDynamics::l1NormConstraintViolation() const {
  return ID_.lpNorm<1>();
}


template <typename MatrixType1, typename MatrixType2>
inline void UnconstrDynamics::getStateFeedbackGain(
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
inline void UnconstrDynamics::getStateFeedbackGain(
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

} // namespace idocp 

#endif // IDOCP_UNCONSTR_DYNAMICS_HXX_ 