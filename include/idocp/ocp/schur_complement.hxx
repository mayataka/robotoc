#ifndef IDOCP_SCHUR_COMPLEMENT_HXX_
#define IDOCP_SCHUR_COMPLEMENT_HXX_

#include "idocp/ocp/schur_complement.hxx"

#include <stdexcept>
#include <assert.h>

#include "Eigen/LU"


namespace idocp {

inline SchurComplement::SchurComplement(const int max_dimA, const int max_dimD)
  : max_dimA_(max_dimA),
    max_dimD_(max_dimD),
    SA_(Eigen::MatrixXd::Zero(max_dimD, max_dimD)),
    SD_(Eigen::MatrixXd::Zero(max_dimA, max_dimA)),
    CAinv_(Eigen::MatrixXd::Zero(max_dimD, max_dimA)),
    BDinv_(Eigen::MatrixXd::Zero(max_dimA, max_dimD)) {
  try {
    if (max_dimA < 0) {
      throw std::out_of_range("invalid value: max_dimA must be non negative!");
    }
    if (max_dimD < 0) {
      throw std::out_of_range("invalid value: max_dimD must be non negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline SchurComplement::SchurComplement() 
  : max_dimA_(0),
    max_dimD_(0),
    SA_(),
    SD_(),
    CAinv_(),
    BDinv_() {
}


inline SchurComplement::~SchurComplement() {
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void SchurComplement::invertWithZeroBottomRightCorner(
    const int dimA, const int dimD, const Eigen::MatrixBase<MatrixType1>& A,
    const Eigen::MatrixBase<MatrixType2>& C,
    const Eigen::MatrixBase<MatrixType3>& Minv) {
  assert(dimA > 0);
  assert(dimA <= max_dimA_);
  assert(dimD > 0);
  assert(dimD <= max_dimD_);
  assert(A.rows() == dimA);
  assert(A.cols() == dimA);
  assert(C.rows() == dimD);
  assert(C.cols() == dimA);
  assert(Minv.rows() == dimA+dimD);
  assert(Minv.cols() == dimA+dimD);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(Minv)
      .topLeftCorner(dimA, dimA).noalias()
      = A.llt().solve(Eigen::MatrixXd::Identity(dimA, dimA));
  CAinv_.topLeftCorner(dimD, dimA).noalias()
      = C * Minv.topLeftCorner(dimA, dimA);
  SA_.topLeftCorner(dimD, dimD).noalias() 
      = CAinv_.topLeftCorner(dimD, dimA) * C.transpose();
  const_cast<Eigen::MatrixBase<MatrixType3>&>(Minv)
      .bottomRightCorner(dimD, dimD).noalias()
      = - SA_.topLeftCorner(dimD, dimD)
             .llt().solve(Eigen::MatrixXd::Identity(dimD, dimD));
  const_cast<Eigen::MatrixBase<MatrixType3>&>(Minv)
      .bottomLeftCorner(dimD, dimA).noalias()
      = - Minv.bottomRightCorner(dimD, dimD) * CAinv_.topLeftCorner(dimD, dimA);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(Minv).topRightCorner(dimA, dimD)
      = Minv.bottomLeftCorner(dimD, dimA).transpose();
  const_cast<Eigen::MatrixBase<MatrixType3>&>(Minv)
      .topLeftCorner(dimA, dimA).noalias()
      -= Minv.bottomLeftCorner(dimD, dimA).transpose()
            * SA_.topLeftCorner(dimD, dimD) * Minv.bottomLeftCorner(dimD, dimA);
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void SchurComplement::invertWithZeroTopLeftCorner(
    const int dimA, const int dimD, const Eigen::MatrixBase<MatrixType1>& B,
    const Eigen::MatrixBase<MatrixType2>& D,
    const Eigen::MatrixBase<MatrixType3>& Minv) {
  assert(dimA > 0);
  assert(dimA <= max_dimA_);
  assert(dimD > 0);
  assert(dimD <= max_dimD_);
  assert(B.rows() == dimA);
  assert(B.cols() == dimD);
  assert(D.rows() == dimD);
  assert(D.cols() == dimD);
  assert(Minv.rows() == dimA+dimD);
  assert(Minv.cols() == dimA+dimD);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(Minv)
      .bottomRightCorner(dimD, dimD).noalias()
      = D.llt().solve(Eigen::MatrixXd::Identity(dimD, dimD));
  BDinv_.topLeftCorner(dimA, dimD).noalias()
      = B * Minv.bottomRightCorner(dimD, dimD);
  SD_.topLeftCorner(dimA, dimA).noalias() 
      = BDinv_.topLeftCorner(dimA, dimD) * B.transpose();
  const_cast<Eigen::MatrixBase<MatrixType3>&>(Minv)
      .topLeftCorner(dimA, dimA).noalias()
      = - SD_.topLeftCorner(dimA, dimA)
             .llt().solve(Eigen::MatrixXd::Identity(dimA, dimA));
  const_cast<Eigen::MatrixBase<MatrixType3>&>(Minv)
      .topRightCorner(dimA, dimD).noalias()
      = - Minv.topLeftCorner(dimA, dimA) * BDinv_.topLeftCorner(dimA, dimD);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(Minv).bottomLeftCorner(dimD, dimA)
      = Minv.topRightCorner(dimA, dimD).transpose();
  const_cast<Eigen::MatrixBase<MatrixType3>&>(Minv)
      .bottomRightCorner(dimD, dimD).noalias()
      -= Minv.topRightCorner(dimA, dimD).transpose()
            * SD_.topLeftCorner(dimA, dimA) * Minv.topRightCorner(dimA, dimD);
}

} // namespace idocp 

#endif // IDOCP_SCHUR_COMPLEMENT_HXX_ 