#ifndef IDOCP_UNCONSTR_KKT_MATRIX_INVERTER_HXX_ 
#define IDOCP_UNCONSTR_KKT_MATRIX_INVERTER_HXX_

#include "idocp/parnmpc/unconstr_kkt_matrix_inverter.hpp"

#include <cassert>


namespace idocp {

inline UnconstrKKTMatrixInverter::UnconstrKKTMatrixInverter(const Robot& robot) 
  : llt_H_(3*robot.dimv()),
    llt_S_(2*robot.dimv()),
    FHinv_(Eigen::MatrixXd::Zero(2*robot.dimv(), 3*robot.dimv())),
    S_(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimH_(3*robot.dimv()), 
    dimkkt_(5*robot.dimv()) {
}


inline UnconstrKKTMatrixInverter::UnconstrKKTMatrixInverter() 
  : llt_H_(),
    llt_S_(),
    FHinv_(),
    S_(),
    dimv_(0), 
    dimx_(0), 
    dimH_(0), 
    dimkkt_(0) {
}


inline UnconstrKKTMatrixInverter::~UnconstrKKTMatrixInverter() {
}


template <typename MatrixType1, typename MatrixType2>
inline void UnconstrKKTMatrixInverter::invert(
    const double dt, const Eigen::MatrixBase<MatrixType1>& H,
    const Eigen::MatrixBase<MatrixType2>& KKT_mat_inv) {
  assert(dt > 0);
  assert(KKT_mat_inv.rows() == dimkkt_);
  assert(KKT_mat_inv.cols() == dimkkt_);
  llt_H_.compute(H);
  assert(llt_H_.info() == Eigen::Success);
  const_cast<Eigen::MatrixBase<MatrixType2>&> (KKT_mat_inv).bottomRightCorner(dimH_, dimH_).noalias()
      = llt_H_.solve(Eigen::MatrixXd::Identity(dimH_, dimH_));
  FHinv_.topRows(dimv_) 
      = - KKT_mat_inv.block(3*dimv_, 2*dimv_, dimv_, 3*dimv_) 
        + dt * KKT_mat_inv.block(4*dimv_, 2*dimv_, dimv_, 3*dimv_);
  FHinv_.bottomRows(dimv_) 
      = dt * KKT_mat_inv.block(2*dimv_, 2*dimv_, dimv_, 3*dimv_) 
        - KKT_mat_inv.block(4*dimv_, 2*dimv_, dimv_, 3*dimv_);
  S_.topLeftCorner(dimv_, dimv_) 
      = - FHinv_.block(0, dimv_, dimv_, dimv_) 
        + dt * FHinv_.block(0, 2*dimv_, dimv_, dimv_);
  S_.topRightCorner(dimv_, dimv_) 
      = dt * FHinv_.block(0, 0, dimv_, dimv_) 
        - FHinv_.block(0, 2*dimv_, dimv_, dimv_);
  S_.bottomLeftCorner(dimv_, dimv_) 
      = - FHinv_.block(dimv_, dimv_, dimv_, dimv_) 
        + dt * FHinv_.block(dimv_, 2*dimv_, dimv_, dimv_);
  S_.bottomRightCorner(dimv_, dimv_) 
      = dt * FHinv_.block(dimv_, 0, dimv_, dimv_) 
        - FHinv_.block(dimv_, 2*dimv_, dimv_, dimv_);
  llt_S_.compute(S_);
  assert(llt_S_.info() == Eigen::Success);
  const_cast<Eigen::MatrixBase<MatrixType2>&> (KKT_mat_inv).topLeftCorner(dimx_, dimx_).noalias()
      = - llt_S_.solve(Eigen::MatrixXd::Identity(dimx_, dimx_));
  const_cast<Eigen::MatrixBase<MatrixType2>&> (KKT_mat_inv).topRightCorner(dimx_, dimH_).noalias()
      = - KKT_mat_inv.topLeftCorner(dimx_, dimx_) * FHinv_;
  const_cast<Eigen::MatrixBase<MatrixType2>&> (KKT_mat_inv).bottomLeftCorner(dimH_, dimx_)
      = KKT_mat_inv.topRightCorner(dimx_, dimH_).transpose();
  const_cast<Eigen::MatrixBase<MatrixType2>&> (KKT_mat_inv).bottomRightCorner(dimH_, dimH_).noalias()
      -= KKT_mat_inv.topRightCorner(dimx_, dimH_).transpose()
          * S_ * KKT_mat_inv.topRightCorner(dimx_, dimH_);
}

} // namespace idocp 

#endif // IDOCP_UNCONSTR_KKT_MATRIX_INVERTER_HXX_ 