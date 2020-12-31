#ifndef IDOCP_SPLIT_UNKKT_MATRIX_HXX_
#define IDOCP_SPLIT_UNKKT_MATRIX_HXX_

#include "idocp/unocp/split_unkkt_matrix.hpp"

#include <cassert>


namespace idocp {

inline SplitUnKKTMatrix::SplitUnKKTMatrix(const Robot& robot) 
  : Q(Eigen::MatrixXd::Zero(3*robot.dimv(), 3*robot.dimv())),
    llt_Q_(3*robot.dimv()),
    llt_S_(2*robot.dimv()),
    FQinv_(Eigen::MatrixXd::Zero(2*robot.dimv(), 3*robot.dimv())),
    S_(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimQ_(3*robot.dimv()), 
    dimKKT_(5*robot.dimv()) {
}


inline SplitUnKKTMatrix::SplitUnKKTMatrix() 
  : Q(),
    llt_Q_(),
    llt_S_(),
    FQinv_(),
    S_(),
    dimv_(0), 
    dimx_(0), 
    dimQ_(0), 
    dimKKT_(0) {
}


inline SplitUnKKTMatrix::~SplitUnKKTMatrix() {
}


inline Eigen::Block<Eigen::MatrixXd> SplitUnKKTMatrix::Qaa() {
  return Q.topLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitUnKKTMatrix::Qaa() const {
  return Q.topLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitUnKKTMatrix::Qaq() {
  return Q.block(0, dimv_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitUnKKTMatrix::Qaq() const {
  return Q.block(0, dimv_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitUnKKTMatrix::Qav() {
  return Q.block(0, dimx_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitUnKKTMatrix::Qav() const {
  return Q.block(0, dimx_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitUnKKTMatrix::Qqa() {
  return Q.block(dimv_, 0, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitUnKKTMatrix::Qqa() const {
  return Q.block(dimv_, 0, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitUnKKTMatrix::Qqq() {
  return Q.block(dimv_, dimv_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitUnKKTMatrix::Qqq() const {
  return Q.block(dimv_, dimv_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitUnKKTMatrix::Qqv() {
  return Q.block(dimv_, dimx_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitUnKKTMatrix::Qqv() const {
  return Q.block(dimv_, dimx_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitUnKKTMatrix::Qva() {
  return Q.block(dimx_, 0, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitUnKKTMatrix::Qva() const {
  return Q.block(dimx_, 0, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitUnKKTMatrix::Qvq() {
  return Q.block(dimx_, dimv_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitUnKKTMatrix::Qvq() const {
  return Q.block(dimx_, dimv_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitUnKKTMatrix::Qvv() {
  return Q.block(dimx_, dimx_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitUnKKTMatrix::Qvv() const {
  return Q.block(dimx_, dimx_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitUnKKTMatrix::Qax() {
  return Q.block(0, dimv_, dimv_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitUnKKTMatrix::Qax() const {
  return Q.block(0, dimv_, dimv_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitUnKKTMatrix::Qxa() {
  return Q.block(dimv_, 0, dimx_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitUnKKTMatrix::Qxa() const {
  return Q.block(dimv_, 0, dimx_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitUnKKTMatrix::Qxx() {
  return Q.block(dimv_, dimv_, dimx_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitUnKKTMatrix::Qxx() const {
  return Q.block(dimv_, dimv_, dimx_, dimx_);
}


inline void SplitUnKKTMatrix::symmetrize() {
  Q.template triangularView<Eigen::StrictlyLower>() 
      = Q.transpose().template triangularView<Eigen::StrictlyLower>();
}


template <typename MatrixType>
inline void SplitUnKKTMatrix::invert(
    const double dtau,
    const Eigen::MatrixBase<MatrixType>& KKT_matrix_inverse) {
  assert(KKT_matrix_inverse.rows() == dimKKT_);
  assert(KKT_matrix_inverse.cols() == dimKKT_);
  llt_Q_.compute(Q);
  assert(llt_Q_.info() == Eigen::Success);
  const_cast<Eigen::MatrixBase<MatrixType>&> (KKT_matrix_inverse).bottomRightCorner(dimQ_, dimQ_).noalias()
      = llt_Q_.solve(Eigen::MatrixXd::Identity(dimQ_, dimQ_));
  FQinv_.topRows(dimv_) 
      = - KKT_matrix_inverse.block(3*dimv_, 2*dimv_, dimv_, 3*dimv_) 
        + dtau * KKT_matrix_inverse.block(4*dimv_, 2*dimv_, dimv_, 3*dimv_);
  FQinv_.bottomRows(dimv_) 
      = dtau * KKT_matrix_inverse.block(2*dimv_, 2*dimv_, dimv_, 3*dimv_) 
        - KKT_matrix_inverse.block(4*dimv_, 2*dimv_, dimv_, 3*dimv_);
  S_.topLeftCorner(dimv_, dimv_) 
      = - FQinv_.block(0, dimv_, dimv_, dimv_) 
        + dtau * FQinv_.block(0, 2*dimv_, dimv_, dimv_);
  S_.topRightCorner(dimv_, dimv_) 
      = dtau * FQinv_.block(0, 0, dimv_, dimv_) 
        - FQinv_.block(0, 2*dimv_, dimv_, dimv_);
  S_.bottomLeftCorner(dimv_, dimv_) 
      = - FQinv_.block(dimv_, dimv_, dimv_, dimv_) 
        + dtau * FQinv_.block(dimv_, 2*dimv_, dimv_, dimv_);
  S_.bottomRightCorner(dimv_, dimv_) 
      = dtau * FQinv_.block(dimv_, 0, dimv_, dimv_) 
        - FQinv_.block(dimv_, 2*dimv_, dimv_, dimv_);
  llt_S_.compute(S_);
  assert(llt_S_.info() == Eigen::Success);
  const_cast<Eigen::MatrixBase<MatrixType>&> (KKT_matrix_inverse).topLeftCorner(dimx_, dimx_).noalias()
      = - llt_S_.solve(Eigen::MatrixXd::Identity(dimx_, dimx_));
  const_cast<Eigen::MatrixBase<MatrixType>&> (KKT_matrix_inverse).topRightCorner(dimx_, dimQ_).noalias()
      = - KKT_matrix_inverse.topLeftCorner(dimx_, dimx_) * FQinv_;
  const_cast<Eigen::MatrixBase<MatrixType>&> (KKT_matrix_inverse).bottomLeftCorner(dimQ_, dimx_)
      = KKT_matrix_inverse.topRightCorner(dimx_, dimQ_).transpose();
  const_cast<Eigen::MatrixBase<MatrixType>&> (KKT_matrix_inverse).bottomRightCorner(dimQ_, dimQ_).noalias()
      -= KKT_matrix_inverse.topRightCorner(dimx_, dimQ_).transpose()
          * S_ * KKT_matrix_inverse.topRightCorner(dimx_, dimQ_);
}


inline void SplitUnKKTMatrix::setZero() {
  Q.setZero();
}


inline int SplitUnKKTMatrix::dimKKT() const {
  return dimKKT_;
}


inline bool SplitUnKKTMatrix::isApprox(const SplitUnKKTMatrix& other) const {
  if (!Q.isApprox(other.Q)) return false;
  else return true;
}


inline bool SplitUnKKTMatrix::hasNaN() const {
  if (Q.hasNaN()) return true;
  else return false;
}

} // namespace idocp 

#endif // IDOCP_SPLIT_UNKKT_MATRIX_HXX_ 