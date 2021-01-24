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
    dimf_(0),
    has_floating_base_(robot.hasFloatingBase()),
    llt_Q_(dimQ_),
    llt_F_(dimx_),
    llt_FPq_(),
    S_full_(Eigen::MatrixXd::Zero(2*robot.dimv()+robot.max_dimf(), 
                                  2*robot.dimv()+robot.max_dimf())),
    Jac_Qinv_full_(Eigen::MatrixXd::Zero(2*robot.dimv()+robot.max_dimf(), 
                                         2*robot.dimv()+robot.dimu())) {
}


inline SplitKKTMatrixInverter::SplitKKTMatrixInverter() 
  : dimv_(0),
    dimu_(0),
    dimx_(0),
    dimQ_(0),
    dimKKT_(0),
    dimf_(0),
    has_floating_base_(false),
    llt_Q_(),
    llt_F_(),
    llt_FPq_(),
    S_full_(),
    Jac_Qinv_full_() {
}


inline SplitKKTMatrixInverter::~SplitKKTMatrixInverter() {
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void SplitKKTMatrixInverter::invert(
    const double dtau, const Eigen::MatrixBase<MatrixType1>& F,
    const Eigen::MatrixBase<MatrixType2>& Q,
    const Eigen::MatrixBase<MatrixType3>& KKT_mat_inv) {
  assert(F.rows() == dimx_);
  assert(F.cols() == dimQ_);
  assert(Q.rows() == dimQ_);
  assert(Q.cols() == dimQ_);
  assert(KKT_mat_inv.rows() == dimKKT_);
  assert(KKT_mat_inv.cols() == dimKKT_);
  llt_Q_.compute(Q);
  assert(llt_Q_.info() == Eigen::Success);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(KKT_mat_inv)
      .bottomRightCorner(dimQ_, dimQ_).noalias()
      = llt_Q_.solve(Eigen::MatrixXd::Identity(dimQ_, dimQ_));
  dimf_ = 0;
  multiplyF(dtau, F, KKT_mat_inv.bottomRightCorner(dimQ_, dimQ_), Jac_Qinv());
  multiplyF(dtau, F, Jac_Qinv().transpose(), S());
  llt_F_.compute(S());
  assert(llt_F_.info() == Eigen::Success);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(KKT_mat_inv)
      .topLeftCorner(dimx_, dimx_).noalias()
      = - llt_F_.solve(Eigen::MatrixXd::Identity(dimx_, dimx_));
  const_cast<Eigen::MatrixBase<MatrixType3>&>(KKT_mat_inv)
      .topRightCorner(dimx_, dimQ_).noalias()
      = - KKT_mat_inv.topLeftCorner(dimx_, dimx_) * Jac_Qinv();
  const_cast<Eigen::MatrixBase<MatrixType3>&>(KKT_mat_inv).bottomLeftCorner(dimQ_, dimx_)
      = KKT_mat_inv.topRightCorner(dimx_, dimQ_).transpose();
  Jac_Qinv().noalias() = S() * KKT_mat_inv.topRightCorner(dimx_, dimQ_);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(KKT_mat_inv)
      .bottomRightCorner(dimQ_, dimQ_).noalias()
      -= KKT_mat_inv.topRightCorner(dimx_, dimQ_).transpose() * Jac_Qinv();
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void SplitKKTMatrixInverter::multiplyF(
    const double dtau, const Eigen::MatrixBase<MatrixType1>& F, 
    const Eigen::MatrixBase<MatrixType2>& mat, 
    const Eigen::MatrixBase<MatrixType3>& res) {
  assert(dtau >= 0);
  assert(F.rows() == dimx_);
  assert(F.cols() == dimQ_);
  if (has_floating_base_) {
    const_cast<Eigen::MatrixBase<MatrixType3>&>(res).template topRows<6>().noalias()
        = F.template block<6, 6>(0, dimu_) * mat.template middleRows<6>(dimu_);
    const_cast<Eigen::MatrixBase<MatrixType3>&>(res).middleRows(6, dimv_-6)
        = - mat.middleRows(dimu_+6, dimv_-6);
    const_cast<Eigen::MatrixBase<MatrixType3>&>(res).template topRows<6>().noalias()
        += F.template block<6, 6>(0, dimu_+dimv_) * mat.template middleRows<6>(dimu_+dimv_);
    const_cast<Eigen::MatrixBase<MatrixType3>&>(res).middleRows(6, dimv_-6).noalias()
        += dtau * mat.middleRows(dimu_+dimv_+6, dimv_-6);
  }
  else {
    const_cast<Eigen::MatrixBase<MatrixType3>&>(res).topRows(dimv_) 
        = - mat.middleRows(dimu_, dimv_);
    const_cast<Eigen::MatrixBase<MatrixType3>&>(res).topRows(dimv_).noalias()
        += dtau * mat.bottomRows(dimv_);
  }
  const_cast<Eigen::MatrixBase<MatrixType3>&>(res).bottomRows(dimv_).noalias()
      = F.bottomRows(dimv_) * mat;
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3, 
          typename MatrixType4>
inline void SplitKKTMatrixInverter::invert(
    const double dtau, const Eigen::MatrixBase<MatrixType1>& F,
    const Eigen::MatrixBase<MatrixType2>& Pq, 
    const Eigen::MatrixBase<MatrixType3>& Q,
    const Eigen::MatrixBase<MatrixType4>& KKT_mat_inv) {
  assert(F.rows() == dimx_);
  assert(F.cols() == dimQ_);
  assert(Pq.cols() == dimv_);
  assert(Q.rows() == dimQ_);
  assert(Q.cols() == dimQ_);
  assert(KKT_mat_inv.rows() == dimKKT_+Pq.rows());
  assert(KKT_mat_inv.cols() == dimKKT_+Pq.rows());
  llt_Q_.compute(Q);
  assert(llt_Q_.info() == Eigen::Success);
  const_cast<Eigen::MatrixBase<MatrixType4>&>(KKT_mat_inv)
      .bottomRightCorner(dimQ_, dimQ_).noalias()
      = llt_Q_.solve(Eigen::MatrixXd::Identity(dimQ_, dimQ_));
  dimf_ = Pq.rows();
  const int dims = dimf_ + dimx_;
  multiplyFPq(dtau, F, Pq, KKT_mat_inv.bottomRightCorner(dimQ_, dimQ_), 
              Jac_Qinv());
  multiplyFPq(dtau, F, Pq, Jac_Qinv().transpose(), S());
  llt_FPq_.compute(S());
  assert(llt_FPq_.info() == Eigen::Success);
  const_cast<Eigen::MatrixBase<MatrixType4>&>(KKT_mat_inv)
      .topLeftCorner(dims, dims).noalias()
      = - llt_FPq_.solve(Eigen::MatrixXd::Identity(dims, dims));
  const_cast<Eigen::MatrixBase<MatrixType4>&>(KKT_mat_inv)
      .topRightCorner(dims, dimQ_).noalias()
      = - KKT_mat_inv.topLeftCorner(dims, dims) * Jac_Qinv();
  const_cast<Eigen::MatrixBase<MatrixType4>&>(KKT_mat_inv).bottomLeftCorner(dimQ_, dims)
      = KKT_mat_inv.topRightCorner(dims, dimQ_).transpose();
  Jac_Qinv().noalias() = S() * KKT_mat_inv.topRightCorner(dims, dimQ_);
  const_cast<Eigen::MatrixBase<MatrixType4>&>(KKT_mat_inv)
      .bottomRightCorner(dimQ_, dimQ_).noalias()
      -= KKT_mat_inv.topRightCorner(dims, dimQ_).transpose() * Jac_Qinv();
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3, 
          typename MatrixType4>
inline void SplitKKTMatrixInverter::multiplyFPq(
    const double dtau, const Eigen::MatrixBase<MatrixType1>& F, 
    const Eigen::MatrixBase<MatrixType2>& Pq, 
    const Eigen::MatrixBase<MatrixType3>& mat, 
    const Eigen::MatrixBase<MatrixType4>& res) {
  assert(dtau >= 0);
  assert(F.rows() == dimx_);
  assert(F.cols() == dimQ_);
  assert(Pq.cols() == dimv_);
  const int dimf = Pq.rows();
  if (has_floating_base_) {
    const_cast<Eigen::MatrixBase<MatrixType4>&>(res).template topRows<6>().noalias()
        = F.template block<6, 6>(0, dimu_) * mat.template middleRows<6>(dimu_);
    const_cast<Eigen::MatrixBase<MatrixType4>&>(res).middleRows(6, dimv_-6)
        = - mat.middleRows(dimu_+6, dimv_-6);
    const_cast<Eigen::MatrixBase<MatrixType4>&>(res).template topRows<6>().noalias()
        += F.template block<6, 6>(0, dimu_+dimv_) * mat.template middleRows<6>(dimu_+dimv_);
    const_cast<Eigen::MatrixBase<MatrixType4>&>(res).middleRows(6, dimv_-6).noalias()
        += dtau * mat.middleRows(dimu_+dimv_+6, dimv_-6);
  }
  else {
    const_cast<Eigen::MatrixBase<MatrixType4>&>(res).topRows(dimv_) 
        = - mat.middleRows(dimu_, dimv_);
    const_cast<Eigen::MatrixBase<MatrixType4>&>(res).topRows(dimv_).noalias()
        += dtau * mat.bottomRows(dimv_);
  }
  const_cast<Eigen::MatrixBase<MatrixType4>&>(res).middleRows(dimv_, dimv_).noalias()
      = F.bottomRows(dimv_) * mat;
  const_cast<Eigen::MatrixBase<MatrixType4>&>(res).bottomRows(dimf).noalias()
      = Pq * mat.middleRows(dimu_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrixInverter::S() {
  return S_full_.topLeftCorner(dimx_+dimf_, dimx_+dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitKKTMatrixInverter::S() const {
  return S_full_.topLeftCorner(dimx_+dimf_, dimx_+dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrixInverter::Jac_Qinv() {
  return Jac_Qinv_full_.topLeftCorner(dimx_+dimf_, dimQ_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitKKTMatrixInverter::Jac_Qinv() const {
  return Jac_Qinv_full_.topLeftCorner(dimx_+dimf_, dimQ_);
}

} // namespace idocp 

#endif // IDOCP_SPLIT_KKT_MATRIX_INVERTER_HXX_ 