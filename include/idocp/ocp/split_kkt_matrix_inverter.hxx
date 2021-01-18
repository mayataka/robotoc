#ifndef IDOCP_SPLIT_KKT_MATRIX_INVERTER_HXX_
#define IDOCP_SPLIT_KKT_MATRIX_INVERTER_HXX_

#include "idocp/ocp/split_kkt_matrix_inverter.hpp"

#include <stdexcept>
#include <assert.h>

#include "Eigen/LU"


namespace idocp {

inline SplitKKTMatrixInverter::SplitKKTMatrixInverter(const Robot& robot)
  : dimv_(robot.dimv()),
    dimu_(robot.dimu()),
    dimx_(2*robot.dimv()),
    dimQ_(2*robot.dimv()+robot.dimu()),
    dimKKT_(4*robot.dimv()+robot.dimu()),
    has_floating_base_(robot.hasFloatingBase()),
    llt_Q_(dimQ_),
    llt_F_(dimx_),
    S_(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    Jac_Qinv_(Eigen::MatrixXd::Zero(2*robot.dimv(), 
                                    2*robot.dimv()+robot.dimu())) {
}


inline SplitKKTMatrixInverter::SplitKKTMatrixInverter() 
  : dimv_(0),
    dimu_(0),
    dimx_(0),
    dimQ_(0),
    dimKKT_(0),
    has_floating_base_(false),
    llt_Q_(),
    llt_F_(),
    S_(),
    Jac_Qinv_() {
}


inline SplitKKTMatrixInverter::~SplitKKTMatrixInverter() {
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void SplitKKTMatrixInverter::invert(
    const double dtau, const Eigen::MatrixBase<MatrixType1>& Jac,
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
  multiplyJac(dtau, Jac, KKT_mat_inv.bottomRightCorner(dimQ_, dimQ_), Jac_Qinv_);
  multiplyJac(dtau, Jac, Jac_Qinv_.transpose(), S_);
  llt_F_.compute(S_);
  assert(llt_F_.info() == Eigen::Success);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(KKT_mat_inv)
      .topLeftCorner(dimx_, dimx_).noalias()
      = - llt_F_.solve(Eigen::MatrixXd::Identity(dimx_, dimx_));
  const_cast<Eigen::MatrixBase<MatrixType3>&>(KKT_mat_inv)
      .topRightCorner(dimx_, dimQ_).noalias()
      = - KKT_mat_inv.topLeftCorner(dimx_, dimx_) * Jac_Qinv_;
  const_cast<Eigen::MatrixBase<MatrixType3>&>(KKT_mat_inv).bottomLeftCorner(dimQ_, dimx_)
      = KKT_mat_inv.topRightCorner(dimx_, dimQ_).transpose();
  Jac_Qinv_.noalias() = S_ * KKT_mat_inv.topRightCorner(dimx_, dimQ_);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(KKT_mat_inv)
      .bottomRightCorner(dimQ_, dimQ_).noalias()
      -= KKT_mat_inv.topRightCorner(dimx_, dimQ_).transpose() * Jac_Qinv_;
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void SplitKKTMatrixInverter::multiplyJac(
    const double dtau, const Eigen::MatrixBase<MatrixType1>& Jac, 
    const Eigen::MatrixBase<MatrixType2>& mat, 
    const Eigen::MatrixBase<MatrixType3>& res) {
  assert(dtau >= 0);
  assert(Jac.rows() == dimx_);
  assert(Jac.cols() == dimQ_);
  if (has_floating_base_) {
    const_cast<Eigen::MatrixBase<MatrixType3>&>(res).topRows(dimv_).noalias()
        = Jac.block(0, dimu_, dimv_, dimv_) * mat.middleRows(dimu_, dimv_);
  }
  else {
    const_cast<Eigen::MatrixBase<MatrixType3>&>(res).topRows(dimv_) 
        = - mat.middleRows(dimu_, dimv_);
  }
  const_cast<Eigen::MatrixBase<MatrixType3>&>(res).topRows(dimv_).noalias()
      += dtau * mat.bottomRows(dimv_);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(res).bottomRows(dimv_).noalias()
      = Jac.bottomRows(dimv_) * mat;
}

} // namespace idocp 

#endif // IDOCP_SPLIT_KKT_MATRIX_INVERTER_HXX_ 