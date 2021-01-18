#ifndef IDOCP_SPLIT_UNKKT_MATRIX_HXX_
#define IDOCP_SPLIT_UNKKT_MATRIX_HXX_

#include "idocp/unocp/split_unkkt_matrix.hpp"

#include <cassert>


namespace idocp {

inline SplitUnKKTMatrix::SplitUnKKTMatrix(const Robot& robot) 
  : Q(Eigen::MatrixXd::Zero(3*robot.dimv(), 3*robot.dimv())),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimKKT_(5*robot.dimv()) {
}


inline SplitUnKKTMatrix::SplitUnKKTMatrix() 
  : Q(),
    dimv_(0), 
    dimx_(0), 
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