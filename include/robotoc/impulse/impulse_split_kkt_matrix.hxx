#ifndef ROBOTOC_IMPULSE_SPLIT_KKT_MATRIX_HXX_
#define ROBOTOC_IMPULSE_SPLIT_KKT_MATRIX_HXX_

#include "robotoc/impulse/impulse_split_kkt_matrix.hpp"

#include <cassert>


namespace robotoc {

inline ImpulseSplitKKTMatrix::ImpulseSplitKKTMatrix(const Robot& robot) 
  : Fxx(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    Qxx(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    Qdvdv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Fqq_prev(),
    Qff_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.max_dimf())),
    Qqf_full_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimi_(0),
    has_floating_base_(robot.hasFloatingBase()) {
  if (robot.hasFloatingBase()) {
    Fqq_prev.resize(robot.dimv(), robot.dimv());
    Fqq_prev.setZero();
  }
}


inline ImpulseSplitKKTMatrix::ImpulseSplitKKTMatrix() 
  : Fxx(),
    Qxx(),
    Qdvdv(),
    Fqq_prev(),
    Qff_full_(),
    Qqf_full_(),
    dimv_(0), 
    dimi_(0) {
}


inline ImpulseSplitKKTMatrix::~ImpulseSplitKKTMatrix() {
}


inline void ImpulseSplitKKTMatrix::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimi_ = impulse_status.dimi();
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Fqq() {
  return Fxx.topLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Fqq() const {
  return Fxx.topLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Fqv() {
  return Fxx.topRightCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Fqv() const {
  return Fxx.topRightCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Fvq() {
  return Fxx.bottomLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Fvq() const {
  return Fxx.bottomLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Fvv() {
  return Fxx.bottomRightCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Fvv() const {
  return Fxx.bottomRightCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Qqq() {
  return Qxx.topLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Qqq() const {
  return Qxx.topLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Qqv() {
  return Qxx.topRightCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Qqv() const {
  return Qxx.topRightCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Qvq() {
  return Qxx.bottomLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Qvq() const {
  return Qxx.bottomLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Qvv() {
  return Qxx.bottomRightCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Qvv() const {
  return Qxx.bottomRightCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Qff() {
  return Qff_full_.topLeftCorner(dimi_, dimi_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Qff() const {
  return Qff_full_.topLeftCorner(dimi_, dimi_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Qqf() {
  return Qqf_full_.topLeftCorner(dimv_, dimi_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Qqf() const {
  return Qqf_full_.topLeftCorner(dimv_, dimi_);
}


inline void ImpulseSplitKKTMatrix::setZero() {
  Fxx.setZero();
  Qxx.setZero();
  Qdvdv.setZero();
  Qff().setZero();
  Qqf().setZero();
  Fqq_prev.setZero();
}


inline int ImpulseSplitKKTMatrix::dimi() const {
  return dimi_;
}

inline bool ImpulseSplitKKTMatrix::isDimensionConsistent() const {
  if (Fxx.cols() != 2*dimv_) return false;
  if (Fxx.rows() != 2*dimv_) return false;
  if (Qxx.cols() != 2*dimv_) return false;
  if (Qxx.rows() != 2*dimv_) return false;
  if (Qdvdv.cols() != dimv_) return false;
  if (Qdvdv.rows() != dimv_) return false;
  if (has_floating_base_) {
    if (Fqq_prev.cols() != dimv_) return false;
    if (Fqq_prev.rows() != dimv_) return false;
  }
  return true;
}


inline bool ImpulseSplitKKTMatrix::isApprox(
    const ImpulseSplitKKTMatrix& other) const {
  if (!Fxx.isApprox(other.Fxx)) return false;
  if (!Qxx.isApprox(other.Qxx)) return false;
  if (!Qdvdv.isApprox(other.Qdvdv)) return false;
  if (!Qff().isApprox(other.Qff())) return false;
  if (!Qqf().isApprox(other.Qqf())) return false;
  if (!Fqq_prev.isApprox(other.Fqq_prev)) return false;
  return true;
}


inline bool ImpulseSplitKKTMatrix::hasNaN() const {
  if (Fxx.hasNaN()) return true;
  if (Qxx.hasNaN()) return true;
  if (Qdvdv.hasNaN()) return true;
  if (Qff().hasNaN()) return true;
  if (Qqf().hasNaN()) return true;
  if (Fqq_prev.hasNaN()) return true;
  return false;
}


inline void ImpulseSplitKKTMatrix::setRandom() {
  Fxx.setRandom();
  const Eigen::MatrixXd Qxx_seed = Eigen::MatrixXd::Random(2*dimv_, 2*dimv_);
  Qxx = Qxx_seed * Qxx_seed.transpose();
  const Eigen::MatrixXd Qdvdvff_seed = Eigen::MatrixXd::Random(dimv_+dimi_, dimv_+dimi_);
  const Eigen::MatrixXd Qdvdvff = Qdvdvff_seed * Qdvdvff_seed.transpose();
  Qdvdv = Qdvdvff.topLeftCorner(dimv_, dimv_);
  Qff() = Qdvdvff.bottomRightCorner(dimi_, dimi_);
  Qqf().setRandom();
  Fqq_prev.setRandom();
}


inline void ImpulseSplitKKTMatrix::setRandom(
    const ImpulseStatus& impulse_status) {
  setImpulseStatus(impulse_status);
  setRandom();
}


inline ImpulseSplitKKTMatrix ImpulseSplitKKTMatrix::Random(
    const Robot& robot) {
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.setRandom();
  return kkt_matrix;
}


inline ImpulseSplitKKTMatrix ImpulseSplitKKTMatrix::Random(
    const Robot& robot, const ImpulseStatus& impulse_status) {
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.setRandom(impulse_status);
  return kkt_matrix;
}

} // namespace robotoc 

#endif // ROBOTOC_IMPULSE_SPLIT_KKT_MATRIX_HXX_ 