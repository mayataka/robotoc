#ifndef IDOCP_SPLIT_KKT_MATRIX_INVERTER_HXX_
#define IDOCP_SPLIT_KKT_MATRIX_INVERTER_HXX_

#include "idocp/ocp/split_kkt_matrix_inverter.hpp"

#include <stdexcept>
#include <assert.h>

#include "Eigen/LU"


namespace idocp {

inline SplitKKTMatrixInverter::SplitKKTMatrixInverter(const Robot& robot)
  : dimx_(2*robot.dimv()),
    dimQ_(2*robot.dimv()+robot.dimu()),
    dimKKT_(4*robot.dimv()+robot.dimu()),
    llt_Q_(dimQ_),
    llt_F_(dimx_),
    S_(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    FQinv_(Eigen::MatrixXd::Zero(2*robot.dimv(), 3*robot.dimv())) {
}


inline SplitKKTMatrixInverter::SplitKKTMatrixInverter() 
  : dimx_(0),
    dimQ_(0),
    dimKKT_(0),
    llt_Q_(),
    llt_F_(),
    S_(),
    FQinv_() {
}


inline SplitKKTMatrixInverter::~SplitKKTMatrixInverter() {
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void SplitKKTMatrixInverter::invert(
    const Eigen::MatrixBase<MatrixType1>& Jac,
    const Eigen::MatrixBase<MatrixType2>& Q,
    const Eigen::MatrixBase<MatrixType3>& KKT_mat_inv) {
  assert(Jac.rows() == dimx_);
  assert(Jac.cols() == dimQ_);
  assert(Q.rows() == dimQ_);
  assert(Q.cols() == dimQ_);
  assert(KKT_mat_inv.rows() == dimKKT_);
  assert(KKT_mat_inv.cols() == dimKKT_);
  llt_Q_.compute(Q);
  assert(llt_Q_.info() == Eigen::Success);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(KKT_mat_inv)
      .bottomRightCorner(dimQ_, dimQ_).noalias()
      = llt_Q_.solve(Eigen::MatrixXd::Identity(dimQ_, dimQ_));
  invert(Jac, const_cast<Eigen::MatrixBase<MatrixType3>&>(KKT_mat_inv));
}


template <typename MatrixType1, typename MatrixType2>
inline void SplitKKTMatrixInverter::invert(
    const Eigen::MatrixBase<MatrixType1>& Jac,
    const Eigen::MatrixBase<MatrixType2>& KKT_mat_inv) {
  assert(Jac.rows() == dimx_);
  assert(Jac.cols() == dimQ_);
  assert(KKT_mat_inv.rows() == dimKKT_);
  assert(KKT_mat_inv.cols() == dimKKT_);
  FQinv_.noalias() = Jac * KKT_mat_inv.bottomRightCorner(dimQ_, dimQ_);
  S_.noalias() = FQinv_ * Jac.transpose();
  llt_F_.compute(S_);
  assert(llt_F_.info() == Eigen::Success);
  const_cast<Eigen::MatrixBase<MatrixType2>&>(KKT_mat_inv)
      .topLeftCorner(dimx_, dimx_).noalias()
      = - llt_F_.solve(Eigen::MatrixXd::Identity(dimx_, dimx_));
  const_cast<Eigen::MatrixBase<MatrixType2>&>(KKT_mat_inv)
      .topRightCorner(dimx_, dimQ_).noalias()
      = - KKT_mat_inv.topLeftCorner(dimx_, dimx_) * FQinv_;
  const_cast<Eigen::MatrixBase<MatrixType2>&>(KKT_mat_inv).bottomLeftCorner(dimQ_, dimx_)
      = KKT_mat_inv.topRightCorner(dimx_, dimQ_).transpose();
  const_cast<Eigen::MatrixBase<MatrixType2>&>(KKT_mat_inv)
      .bottomRightCorner(dimQ_, dimQ_).noalias()
      -= KKT_mat_inv.topRightCorner(dimx_, dimQ_).transpose()
            * S_ * KKT_mat_inv.topRightCorner(dimx_, dimQ_);
}

} // namespace idocp 

#endif // IDOCP_SPLIT_KKT_MATRIX_INVERTER_HXX_ 