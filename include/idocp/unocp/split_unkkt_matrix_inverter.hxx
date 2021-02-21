#ifndef IDOCP_SPLIT_UNKKT_MATRIX_INVERTER_HXX_
#define IDOCP_SPLIT_UNKKT_MATRIX_INVERTER_HXX_

#include "idocp/unocp/split_unkkt_matrix_inverter.hpp"

#include <cassert>


namespace idocp {

inline SplitUnKKTMatrixInverter::SplitUnKKTMatrixInverter(const Robot& robot) 
  : llt_Q_(3*robot.dimv()),
    llt_S_(2*robot.dimv()),
    FQinv_(Eigen::MatrixXd::Zero(2*robot.dimv(), 3*robot.dimv())),
    S_(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimQ_(3*robot.dimv()), 
    dimKKT_(5*robot.dimv()) {
}


inline SplitUnKKTMatrixInverter::SplitUnKKTMatrixInverter() 
  : llt_Q_(),
    llt_S_(),
    FQinv_(),
    S_(),
    dimv_(0), 
    dimx_(0), 
    dimQ_(0), 
    dimKKT_(0) {
}


inline SplitUnKKTMatrixInverter::~SplitUnKKTMatrixInverter() {
}


template <typename MatrixType1, typename MatrixType2>
inline void SplitUnKKTMatrixInverter::invert(
    const double dt, const Eigen::MatrixBase<MatrixType1>& Q,
    const Eigen::MatrixBase<MatrixType2>& KKT_matrix_inverse) {
  assert(dt > 0);
  assert(KKT_matrix_inverse.rows() == dimKKT_);
  assert(KKT_matrix_inverse.cols() == dimKKT_);
  llt_Q_.compute(Q);
  assert(llt_Q_.info() == Eigen::Success);
  const_cast<Eigen::MatrixBase<MatrixType2>&> (KKT_matrix_inverse).bottomRightCorner(dimQ_, dimQ_).noalias()
      = llt_Q_.solve(Eigen::MatrixXd::Identity(dimQ_, dimQ_));
  FQinv_.topRows(dimv_) 
      = - KKT_matrix_inverse.block(3*dimv_, 2*dimv_, dimv_, 3*dimv_) 
        + dt * KKT_matrix_inverse.block(4*dimv_, 2*dimv_, dimv_, 3*dimv_);
  FQinv_.bottomRows(dimv_) 
      = dt * KKT_matrix_inverse.block(2*dimv_, 2*dimv_, dimv_, 3*dimv_) 
        - KKT_matrix_inverse.block(4*dimv_, 2*dimv_, dimv_, 3*dimv_);
  S_.topLeftCorner(dimv_, dimv_) 
      = - FQinv_.block(0, dimv_, dimv_, dimv_) 
        + dt * FQinv_.block(0, 2*dimv_, dimv_, dimv_);
  S_.topRightCorner(dimv_, dimv_) 
      = dt * FQinv_.block(0, 0, dimv_, dimv_) 
        - FQinv_.block(0, 2*dimv_, dimv_, dimv_);
  S_.bottomLeftCorner(dimv_, dimv_) 
      = - FQinv_.block(dimv_, dimv_, dimv_, dimv_) 
        + dt * FQinv_.block(dimv_, 2*dimv_, dimv_, dimv_);
  S_.bottomRightCorner(dimv_, dimv_) 
      = dt * FQinv_.block(dimv_, 0, dimv_, dimv_) 
        - FQinv_.block(dimv_, 2*dimv_, dimv_, dimv_);
  llt_S_.compute(S_);
  assert(llt_S_.info() == Eigen::Success);
  const_cast<Eigen::MatrixBase<MatrixType2>&> (KKT_matrix_inverse).topLeftCorner(dimx_, dimx_).noalias()
      = - llt_S_.solve(Eigen::MatrixXd::Identity(dimx_, dimx_));
  const_cast<Eigen::MatrixBase<MatrixType2>&> (KKT_matrix_inverse).topRightCorner(dimx_, dimQ_).noalias()
      = - KKT_matrix_inverse.topLeftCorner(dimx_, dimx_) * FQinv_;
  const_cast<Eigen::MatrixBase<MatrixType2>&> (KKT_matrix_inverse).bottomLeftCorner(dimQ_, dimx_)
      = KKT_matrix_inverse.topRightCorner(dimx_, dimQ_).transpose();
  const_cast<Eigen::MatrixBase<MatrixType2>&> (KKT_matrix_inverse).bottomRightCorner(dimQ_, dimQ_).noalias()
      -= KKT_matrix_inverse.topRightCorner(dimx_, dimQ_).transpose()
          * S_ * KKT_matrix_inverse.topRightCorner(dimx_, dimQ_);
}

} // namespace idocp 

#endif // IDOCP_SPLIT_UNKKT_MATRIX_INVERTER_HXX_ 