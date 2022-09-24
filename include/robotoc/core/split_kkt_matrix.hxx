#ifndef ROBOTOC_SPLIT_KKT_MATRIX_HXX_
#define ROBOTOC_SPLIT_KKT_MATRIX_HXX_

#include "robotoc/core/split_kkt_matrix.hpp"

#include <cassert>


namespace robotoc {

inline void SplitKKTMatrix::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
}


inline void SplitKKTMatrix::setContactStatus(
    const ImpulseStatus& contact_status) {
  dimf_ = contact_status.dimf();
}


inline void SplitKKTMatrix::setSwitchingConstraintDimension(const int dims) {
  assert(dims >= 0);
  assert(dims <= Phit_full_.size());
  dims_ = dims;
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Fqq() {
  return Fxx.topLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Fqq() const {
  return Fxx.topLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Fqv() {
  return Fxx.topRightCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Fqv() const {
  return Fxx.topRightCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Fvq() {
  return Fxx.bottomLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Fvq() const {
  return Fxx.bottomLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Fvv() {
  return Fxx.bottomRightCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Fvv() const {
  return Fxx.bottomRightCorner(dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTMatrix::fq() {
  return fx.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> SplitKKTMatrix::fq() const {
  return fx.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTMatrix::fv() {
  return fx.tail(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> SplitKKTMatrix::fv() const {
  return fx.tail(dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Phix() {
  return Phix_full_.topLeftCorner(dims_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Phix() const {
  return Phix_full_.topLeftCorner(dims_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Phiq() {
  return Phix_full_.topLeftCorner(dims_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Phiq() const {
  return Phix_full_.topLeftCorner(dims_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Phiv() {
  return Phix_full_.topRightCorner(dims_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Phiv() const {
  return Phix_full_.topRightCorner(dims_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Phia() {
  return Phia_full_.topLeftCorner(dims_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Phia() const {
  return Phia_full_.topLeftCorner(dims_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Phiu() {
  return Phiu_full_.topLeftCorner(dims_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Phiu() const {
  return Phiu_full_.topLeftCorner(dims_, dimu_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTMatrix::Phit() {
  return Phit_full_.head(dims_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> SplitKKTMatrix::Phit() const {
  return Phit_full_.head(dims_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qqq() {
  return Qxx.topLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qqq() const {
  return Qxx.topLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qqv() {
  return Qxx.topRightCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qqv() const {
  return Qxx.topRightCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qvq() {
  return Qxx.bottomLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qvq() const {
  return Qxx.bottomLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qvv() {
  return Qxx.bottomRightCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qvv() const {
  return Qxx.bottomRightCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qqu() {
  return Qxu.topLeftCorner(dimv_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qqu() const {
  return Qxu.topLeftCorner(dimv_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qvu() {
  return Qxu.bottomLeftCorner(dimv_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qvu() const {
  return Qxu.bottomLeftCorner(dimv_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qff() {
  return Qff_full_.topLeftCorner(dimf_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qff() const {
  return Qff_full_.topLeftCorner(dimf_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qqf() {
  return Qqf_full_.topLeftCorner(dimv_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qqf() const {
  return Qqf_full_.topLeftCorner(dimv_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTMatrix::hq() {
  return hx.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> SplitKKTMatrix::hq() const {
  return hx.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTMatrix::hv() {
  return hx.tail(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> SplitKKTMatrix::hv() const {
  return hx.tail(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTMatrix::hf() {
  return hf_full_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> SplitKKTMatrix::hf() const {
  return hf_full_.head(dimf_);
}


inline void SplitKKTMatrix::setZero() {
  Fxx.setZero();
  Fvu.setZero();
  fx.setZero();
  Phix().setZero();
  Phia().setZero();
  Phiu().setZero();
  Phit().setZero();
  Qxx.setZero();
  Qaa.setZero();
  Qdvdv.setZero();
  Qxu.setZero();
  Quu.setZero();
  Qff().setZero();
  Qqf().setZero();
  Qtt = 0;
  Qtt_prev = 0;
  hx.setZero();
  hu.setZero();
  ha.setZero();
  hf().setZero();
}


inline int SplitKKTMatrix::dimf() const {
  return dimf_;
}


inline int SplitKKTMatrix::dims() const {
  return dims_;
}

} // namespace robotoc 

#endif // ROBOTOC_SPLIT_KKT_MATRIX_HXX_ 