#ifndef IDOCP_IMPULSE_SPLIT_KKT_MATRIX_INVERTER_HXX_
#define IDOCP_IMPULSE_SPLIT_KKT_MATRIX_INVERTER_HXX_

#include "idocp/impulse/impulse_split_kkt_matrix_inverter.hpp"

#include <iostream>
#include <stdexcept>
#include <cassert>


namespace idocp {

inline ImpulseSplitKKTMatrixInverter::ImpulseSplitKKTMatrixInverter(
    const int dimv, const int max_dimf)
  : dimv_(dimv),
    max_dimf_(max_dimf),
    S_(Eigen::MatrixXd::Zero(2*dimv+max_dimf_, 2*dimv+max_dimf_)),
    Jac_Qinv_(Eigen::MatrixXd::Zero(2*dimv+max_dimf_, 2*dimv+max_dimf_)),
    llt_() {
  try {
    if (dimv <= 0) {
      throw std::out_of_range("invalid value: dimv must be non negative!");
    }
    if (max_dimf <= 0) {
      throw std::out_of_range("invalid value: max_dimf must be non negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline ImpulseSplitKKTMatrixInverter::ImpulseSplitKKTMatrixInverter() 
  : dimv_(0),
    max_dimf_(0),
    S_(),
    Jac_Qinv_(),
    llt_() {
}


inline ImpulseSplitKKTMatrixInverter::~ImpulseSplitKKTMatrixInverter() {
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void ImpulseSplitKKTMatrixInverter::invert(
    const Eigen::MatrixBase<MatrixType1>& Jac,
    const Eigen::MatrixBase<MatrixType2>& Q,
    const Eigen::MatrixBase<MatrixType3>& KKT_mat_inv) {
  const int dimQ = Jac.rows();
  assert(dimQ <= 2*dimv_+max_dimf_);
  assert(dimQ <= 2*dimv_+max_dimf_);
  assert(Jac.rows() == dimQ);
  assert(Jac.cols() == dimQ);
  assert(Q.rows() == dimQ);
  assert(Q.cols() == dimQ);
  assert(KKT_mat_inv.rows() == 2*dimQ);
  assert(KKT_mat_inv.cols() == 2*dimQ);
  llt_.compute(Q);
  assert(llt_.info() == Eigen::Success);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(KKT_mat_inv)
      .bottomRightCorner(dimQ, dimQ).noalias()
      = llt_.solve(Eigen::MatrixXd::Identity(dimQ, dimQ));
  multiplyJac(Jac, KKT_mat_inv.bottomRightCorner(dimQ, dimQ), 
             Jac_Qinv_.topLeftCorner(dimQ, dimQ));
  multiplyJac(Jac, Jac_Qinv_.topLeftCorner(dimQ, dimQ).transpose(),
             S_.topLeftCorner(dimQ, dimQ));
  llt_.compute(S_.topLeftCorner(dimQ, dimQ));
  assert(llt_.info() == Eigen::Success);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(KKT_mat_inv)
      .topLeftCorner(dimQ, dimQ).noalias()
      = - llt_.solve(Eigen::MatrixXd::Identity(dimQ, dimQ));
  const_cast<Eigen::MatrixBase<MatrixType3>&>(KKT_mat_inv)
      .topRightCorner(dimQ, dimQ).noalias()
      = - KKT_mat_inv.topLeftCorner(dimQ, dimQ) 
            * Jac_Qinv_.topLeftCorner(dimQ, dimQ);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(KKT_mat_inv).bottomLeftCorner(dimQ, dimQ)
      = KKT_mat_inv.topRightCorner(dimQ, dimQ).transpose();
  Jac_Qinv_.topLeftCorner(dimQ, dimQ).noalias()
      = S_.topLeftCorner(dimQ, dimQ) * KKT_mat_inv.topRightCorner(dimQ, dimQ);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(KKT_mat_inv)
      .bottomRightCorner(dimQ, dimQ).noalias()
      -= KKT_mat_inv.topRightCorner(dimQ, dimQ).transpose()
            * Jac_Qinv_.topLeftCorner(dimQ, dimQ);
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
    const_cast<Eigen::MatrixBase<MatrixType3>&>(res).topRows(dimv_).noalias()
        = Jac.block(0, dimf, dimv_, dimv_) * mat.middleRows(dimf, dimv_);
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

} // namespace idocp 

#endif // IDOCP_IMPULSE_SPLIT_KKT_MATRIX_INVERTER_HXX_