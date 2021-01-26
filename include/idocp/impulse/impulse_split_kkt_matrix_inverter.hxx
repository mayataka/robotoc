#ifndef IDOCP_IMPULSE_SPLIT_KKT_MATRIX_INVERTER_HXX_
#define IDOCP_IMPULSE_SPLIT_KKT_MATRIX_INVERTER_HXX_

#include "idocp/impulse/impulse_split_kkt_matrix_inverter.hpp"

#include <iostream>
#include <cassert>


namespace idocp {

inline ImpulseSplitKKTMatrixInverter::ImpulseSplitKKTMatrixInverter(
    const Robot& robot)
  : dimv_(robot.dimv()),
    max_dimf_(robot.max_dimf()),
    dimQ_(2*robot.dimv()),
    has_floating_base_(robot.hasFloatingBase()),
    regularization_(false),
    reg_(1.0e-09),
    S_full_(Eigen::MatrixXd::Zero(2*robot.dimv()+robot.max_dimf(), 
                                  2*robot.dimv()+robot.max_dimf())),
    Jac_Qinv_full_(Eigen::MatrixXd::Zero(2*robot.dimv()+robot.max_dimf(), 
                                         2*robot.dimv()+robot.max_dimf())),
    llt_() {
}


inline ImpulseSplitKKTMatrixInverter::ImpulseSplitKKTMatrixInverter() 
  : dimv_(0),
    max_dimf_(0),
    dimQ_(0),
    has_floating_base_(false),
    regularization_(false),
    reg_(1.0e-09),
    S_full_(),
    Jac_Qinv_full_(),
    llt_() {
}


inline ImpulseSplitKKTMatrixInverter::~ImpulseSplitKKTMatrixInverter() {
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void ImpulseSplitKKTMatrixInverter::invert(
    const Eigen::MatrixBase<MatrixType1>& Jac,
    const Eigen::MatrixBase<MatrixType2>& Q,
    const Eigen::MatrixBase<MatrixType3>& KKT_mat_inv) {
  dimQ_ = Jac.rows();
  assert(dimQ_ <= 2*dimv_+max_dimf_);
  assert(dimQ_ <= 2*dimv_+max_dimf_);
  assert(Jac.rows() == dimQ_);
  assert(Jac.cols() == dimQ_);
  assert(Q.rows() == dimQ_);
  assert(Q.cols() == dimQ_);
  assert(KKT_mat_inv.rows() == 2*dimQ_);
  assert(KKT_mat_inv.cols() == 2*dimQ_);
  llt_ = Eigen::LLT<Eigen::MatrixXd>(dimQ_);
  llt_.compute(Q);
  assert(llt_.info() == Eigen::Success);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(KKT_mat_inv)
      .bottomRightCorner(dimQ_, dimQ_).noalias()
      = llt_.solve(Eigen::MatrixXd::Identity(dimQ_, dimQ_));
  multiplyJac(Jac, KKT_mat_inv.bottomRightCorner(dimQ_, dimQ_), Jac_Qinv());
  multiplyJac(Jac, Jac_Qinv().transpose(), S());
  if (regularization_) {
    const int dimf = dimQ_ - 2*dimv_;
    S().diagonal().tail(dimf).array() += reg_;
  }
  llt_.compute(S());
  assert(llt_.info() == Eigen::Success);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(KKT_mat_inv)
      .topLeftCorner(dimQ_, dimQ_).noalias()
      = - llt_.solve(Eigen::MatrixXd::Identity(dimQ_, dimQ_));
  const_cast<Eigen::MatrixBase<MatrixType3>&>(KKT_mat_inv)
      .topRightCorner(dimQ_, dimQ_).noalias()
      = - KKT_mat_inv.topLeftCorner(dimQ_, dimQ_) * Jac_Qinv();
  const_cast<Eigen::MatrixBase<MatrixType3>&>(KKT_mat_inv).bottomLeftCorner(dimQ_, dimQ_)
      = KKT_mat_inv.topRightCorner(dimQ_, dimQ_).transpose();
  Jac_Qinv().noalias() = S() * KKT_mat_inv.topRightCorner(dimQ_, dimQ_);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(KKT_mat_inv)
      .bottomRightCorner(dimQ_, dimQ_).noalias()
      -= KKT_mat_inv.topRightCorner(dimQ_, dimQ_).transpose() * Jac_Qinv();
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void ImpulseSplitKKTMatrixInverter::multiplyJac(
    const Eigen::MatrixBase<MatrixType1>& Jac, 
    const Eigen::MatrixBase<MatrixType2>& mat, 
    const Eigen::MatrixBase<MatrixType3>& res) {
  const int dimQ = Jac.rows();
  const int dimf = dimQ - 2*dimv_;
  assert(dimQ <= 2*dimv_+max_dimf_);
  assert(dimQ <= 2*dimv_+max_dimf_);
  assert(Jac.rows() == dimQ);
  assert(Jac.cols() == dimQ);
  assert(mat.rows() == dimQ);
  assert(mat.cols() == dimQ);
  assert(res.rows() == dimQ);
  assert(res.cols() == dimQ);
  if (has_floating_base_) {
    const_cast<Eigen::MatrixBase<MatrixType3>&>(res).template topRows<6>().noalias()
        = Jac.template block<6, 6>(0, dimf) * mat.template middleRows<6>(dimf);
    const_cast<Eigen::MatrixBase<MatrixType3>&>(res).middleRows(6, dimv_-6)
        = - mat.middleRows(dimf+6, dimv_-6);
  }
  else {
    const_cast<Eigen::MatrixBase<MatrixType3>&>(res).topRows(dimv_) 
        = - mat.middleRows(dimf, dimv_);
  }
  const_cast<Eigen::MatrixBase<MatrixType3>&>(res).middleRows(dimv_, dimv_).noalias()
      = Jac.block(dimv_, 0, dimv_, dimv_+dimf) * mat.topRows(dimv_+dimf);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(res).middleRows(dimv_, dimv_).noalias() 
      -= mat.bottomRows(dimv_);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(res).bottomRows(dimf).noalias()
      = Jac.block(2*dimv_, dimf, dimf, 2*dimv_) * mat.bottomRows(2*dimv_);
}


inline void ImpulseSplitKKTMatrixInverter::enableRegularization(
    const double reg) {
  assert(reg >= 0);
  regularization_ = true;
  reg_ = reg;
}


inline void ImpulseSplitKKTMatrixInverter::disableRegularization() {
  regularization_ = false;
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrixInverter::S() {
  return S_full_.topLeftCorner(dimQ_, dimQ_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrixInverter::S() const {
  return S_full_.topLeftCorner(dimQ_, dimQ_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrixInverter::Jac_Qinv() {
  return Jac_Qinv_full_.topLeftCorner(dimQ_, dimQ_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrixInverter::Jac_Qinv() const {
  return Jac_Qinv_full_.topLeftCorner(dimQ_, dimQ_);
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_SPLIT_KKT_MATRIX_INVERTER_HXX_