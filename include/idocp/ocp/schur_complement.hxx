#ifndef IDOCP_SCHUR_COMPLEMENT_HXX_
#define IDOCP_SCHUR_COMPLEMENT_HXX_

#include "idocp/ocp/schur_complement.hxx"

#include <stdexcept>
#include <assert.h>

#include "Eigen/LU"


namespace idocp {

inline SchurComplement::SchurComplement(const int dimA, const int dimD)
  : dimA_(dimA),
    dimD_(dimD),
    llt_Ainv_(dimA),
    llt_Dinv_(dimD),
    SA_(Eigen::MatrixXd::Zero(dimD, dimD)),
    SD_(Eigen::MatrixXd::Zero(dimA, dimA)),
    CAinv_(Eigen::MatrixXd::Zero(dimD, dimA)),
    BDinv_(Eigen::MatrixXd::Zero(dimA, dimD)) {
  try {
    if (dimA < 0) {
      throw std::out_of_range("invalid value: dimA must be non negative!");
    }
    if (dimD < 0) {
      throw std::out_of_range("invalid value: dimD must be non negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline SchurComplement::SchurComplement() 
  : dimA_(0),
    dimD_(0),
    SA_(),
    SD_(),
    CAinv_(),
    BDinv_() {
}


inline SchurComplement::~SchurComplement() {
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void SchurComplement::invertWithZeroBottomRightCorner(
    const Eigen::MatrixBase<MatrixType1>& A,
    const Eigen::MatrixBase<MatrixType2>& C,
    const Eigen::MatrixBase<MatrixType3>& Minv) {
  assert(A.rows() == dimA_);
  assert(A.cols() == dimA_);
  assert(C.rows() == dimD_);
  assert(C.cols() == dimA_);
  assert(Minv.rows() == dimA_+dimD_);
  assert(Minv.cols() == dimA_+dimD_);
  llt_Ainv_.compute(A);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(Minv)
      .topLeftCorner(dimA_, dimA_).noalias()
      = llt_Ainv_.solve(Eigen::MatrixXd::Identity(dimA_, dimA_));
  invertWithZeroBottomRightCorner(
      C, const_cast<Eigen::MatrixBase<MatrixType3>&>(Minv));
}


template <typename MatrixType1, typename MatrixType2>
inline void SchurComplement::invertWithZeroBottomRightCorner(
    const Eigen::MatrixBase<MatrixType1>& C,
    const Eigen::MatrixBase<MatrixType2>& Minv) {
  assert(C.rows() == dimD_);
  assert(C.cols() == dimA_);
  assert(Minv.rows() == dimA_+dimD_);
  assert(Minv.cols() == dimA_+dimD_);
  CAinv_.noalias() = C * Minv.topLeftCorner(dimA_, dimA_);
  SA_.noalias() = CAinv_ * C.transpose();
  llt_Dinv_.compute(SA_);
  const_cast<Eigen::MatrixBase<MatrixType2>&>(Minv)
      .bottomRightCorner(dimD_, dimD_).noalias()
      = - llt_Dinv_.solve(Eigen::MatrixXd::Identity(dimD_, dimD_));
  const_cast<Eigen::MatrixBase<MatrixType2>&>(Minv)
      .bottomLeftCorner(dimD_, dimA_).noalias()
      = - Minv.bottomRightCorner(dimD_, dimD_) * CAinv_;
  const_cast<Eigen::MatrixBase<MatrixType2>&>(Minv).topRightCorner(dimA_, dimD_)
      = Minv.bottomLeftCorner(dimD_, dimA_).transpose();
  const_cast<Eigen::MatrixBase<MatrixType2>&>(Minv)
      .topLeftCorner(dimA_, dimA_).noalias()
      -= Minv.bottomLeftCorner(dimD_, dimA_).transpose()
            * SA_ * Minv.bottomLeftCorner(dimD_, dimA_);
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void SchurComplement::invertWithZeroTopLeftCorner(
    const Eigen::MatrixBase<MatrixType1>& B,
    const Eigen::MatrixBase<MatrixType2>& D,
    const Eigen::MatrixBase<MatrixType3>& Minv) {
  assert(B.rows() == dimA_);
  assert(B.cols() == dimD_);
  assert(D.rows() == dimD_);
  assert(D.cols() == dimD_);
  assert(Minv.rows() == dimA_+dimD_);
  assert(Minv.cols() == dimA_+dimD_);
  llt_Dinv_.compute(D);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(Minv)
      .bottomRightCorner(dimD_, dimD_).noalias()
      = llt_Dinv_.solve(Eigen::MatrixXd::Identity(dimD_, dimD_));
  invertWithZeroTopLeftCorner( 
      B, const_cast<Eigen::MatrixBase<MatrixType3>&>(Minv));
}


template <typename MatrixType1, typename MatrixType2>
inline void SchurComplement::invertWithZeroTopLeftCorner(
    const Eigen::MatrixBase<MatrixType1>& B,
    const Eigen::MatrixBase<MatrixType2>& Minv) {
  assert(B.rows() == dimA_);
  assert(B.cols() == dimD_);
  assert(Minv.rows() == dimA_+dimD_);
  assert(Minv.cols() == dimA_+dimD_);
  BDinv_.noalias() = B * Minv.bottomRightCorner(dimD_, dimD_);
  SD_.noalias() = BDinv_ * B.transpose();
  llt_Ainv_.compute(SD_);
  const_cast<Eigen::MatrixBase<MatrixType2>&>(Minv)
      .topLeftCorner(dimA_, dimA_).noalias()
      = - llt_Ainv_.solve(Eigen::MatrixXd::Identity(dimA_, dimA_));
  const_cast<Eigen::MatrixBase<MatrixType2>&>(Minv)
      .topRightCorner(dimA_, dimD_).noalias()
      = - Minv.topLeftCorner(dimA_, dimA_) * BDinv_;
  const_cast<Eigen::MatrixBase<MatrixType2>&>(Minv).bottomLeftCorner(dimD_, dimA_)
      = Minv.topRightCorner(dimA_, dimD_).transpose();
  const_cast<Eigen::MatrixBase<MatrixType2>&>(Minv)
      .bottomRightCorner(dimD_, dimD_).noalias()
      -= Minv.topRightCorner(dimA_, dimD_).transpose()
            * SD_ * Minv.topRightCorner(dimA_, dimD_);
}

} // namespace idocp 

#endif // IDOCP_SCHUR_COMPLEMENT_HXX_ 