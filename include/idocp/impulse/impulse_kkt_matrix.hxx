#ifndef IDOCP_IMPULSE_KKT_MATRIX_HXX_
#define IDOCP_IMPULSE_KKT_MATRIX_HXX_

#include "idocp/impulse/impulse_kkt_matrix.hpp"

#include <cassert>


namespace idocp {

inline ImpulseKKTMatrix::ImpulseKKTMatrix(const Robot& robot) 
  : Fqq_prev(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    schur_complement_(2*robot.dimv(), 2*robot.dimv()+robot.dimu()),
    F_(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv()+robot.max_dimf())),
    C_(Eigen::MatrixXd::Zero(2*robot.max_dimf(), 
                             2*robot.dimv()+robot.max_dimf())),
    Q_(Eigen::MatrixXd::Zero(3*robot.dimv()+robot.max_dimf(), 
                             3*robot.dimv()+robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimf_(0), 
    q_begin_(robot.dimv()),
    v_begin_(2*robot.dimv()),
    dimKKT_(2*robot.dimv()) {
}


inline ImpulseKKTMatrix::ImpulseKKTMatrix() 
  : Fqq_prev(),
    schur_complement_(),
    F_(),
    C_(),
    Q_(),
    dimv_(0), 
    dimx_(0), 
    dimf_(0), 
    q_begin_(0),
    v_begin_(0),
    dimKKT_(0) {
}


inline ImpulseKKTMatrix::~ImpulseKKTMatrix() {
}


inline void ImpulseKKTMatrix::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimf_ = impulse_status.dimp();
  q_begin_ = dimv_ + dimf_;
  v_begin_ = 2*dimv_ + dimf_;
  dimKKT_ = 2*dimv_ + dimf_;
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Fqf() {
  return F_.block(0, 0, dimv_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> ImpulseKKTMatrix::Fqf() const {
  return F_.block(0, 0, dimv_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Fqq() {
  return F_.block(0, dimf_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> ImpulseKKTMatrix::Fqq() const {
  return F_.block(0, dimf_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Fqv() {
  return F_.block(0, dimf_+dimv_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> ImpulseKKTMatrix::Fqv() const {
  return F_.block(0, dimf_+dimv_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Fvf() {
  return F_.block(dimv_, 0, dimv_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> ImpulseKKTMatrix::Fvf() const {
  return F_.block(dimv_, 0, dimv_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Fvq() {
  return F_.block(dimv_, dimf_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> ImpulseKKTMatrix::Fvq() const {
  return F_.block(dimv_, dimf_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Fvv() {
  return F_.block(dimv_, dimf_+dimv_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> ImpulseKKTMatrix::Fvv() const {
  return F_.block(dimv_, dimf_+dimv_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Fxf() {
  return F_.block(0, 0, dimx_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> ImpulseKKTMatrix::Fxf() const {
  return F_.block(0, 0, dimx_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Fxx() {
  return F_.block(0, dimf_, dimx_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> ImpulseKKTMatrix::Fxx() const {
  return F_.block(0, dimf_, dimx_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Pq() {
  return C_.block(0, 0, dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> ImpulseKKTMatrix::Pq() const {
  return C_.block(0, 0, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Vq() {
  return C_.block(dimf_, 0, dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> ImpulseKKTMatrix::Vq() const {
  return C_.block(dimf_, 0, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Vv() {
  return C_.block(dimf_, dimv_, dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> ImpulseKKTMatrix::Vv() const {
  return C_.block(dimf_, dimv_, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qdvdvff() {
  return Q_.topLeftCorner(dimv_+dimf_, dimv_+dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseKKTMatrix::Qdvdvff() const {
  return Q_.topLeftCorner(dimv_+dimf_, dimv_+dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qdvdv() {
  return Q_.topLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseKKTMatrix::Qdvdv() const {
  return Q_.topLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qff() {
  return Q_.block(dimv_, dimv_, dimf_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> ImpulseKKTMatrix::Qff() const {
  return Q_.block(dimv_, dimv_, dimf_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qfq() {
  return Q_.block(dimv_, dimv_+dimf_, dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> ImpulseKKTMatrix::Qfq() const {
  return Q_.block(dimv_, dimv_+dimf_, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qqf() {
  return Q_.block(dimv_+dimf_, dimv_, dimv_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> ImpulseKKTMatrix::Qqf() const {
  return Q_.block(dimv_+dimf_, dimv_, dimv_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qqq() {
  return Q_.block(q_begin_, q_begin_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> ImpulseKKTMatrix::Qqq() const {
  return Q_.block(q_begin_, q_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qqv() {
  return Q_.block(q_begin_, v_begin_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> ImpulseKKTMatrix::Qqv() const {
  return Q_.block(q_begin_, v_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qvq() {
  return Q_.block(v_begin_, q_begin_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> ImpulseKKTMatrix::Qvq() const {
  return Q_.block(v_begin_, q_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qvv() {
  return Q_.block(v_begin_, v_begin_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> ImpulseKKTMatrix::Qvv() const {
  return Q_.block(v_begin_, v_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qxx() {
  return Q_.block(q_begin_, q_begin_, dimx_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> ImpulseKKTMatrix::Qxx() const {
  return Q_.block(q_begin_, q_begin_, dimx_, dimx_);
}


inline void ImpulseKKTMatrix::symmetrize() {
  Q_.template triangularView<Eigen::StrictlyLower>() 
      = Q_.transpose().template triangularView<Eigen::StrictlyLower>();
}


template <typename MatrixType>
inline void ImpulseKKTMatrix::invert(
    const Eigen::MatrixBase<MatrixType>& KKT_matrix_inverse) {
}


inline void ImpulseKKTMatrix::setZero() {
  Fqq_prev.setZero();
  F_.setZero();
  C_.setZero();
  Q_.setZero();
}


inline int ImpulseKKTMatrix::dimKKT() const {
  return dimKKT_;
}


inline int ImpulseKKTMatrix::dimf() const {
  return dimf_;
}


inline bool ImpulseKKTMatrix::isApprox(const ImpulseKKTMatrix& other) const {
  if (!Fxf().isApprox(other.Fxf())) return false;
  if (!Fxx().isApprox(other.Fxx())) return false;
  if (!Pq().isApprox(other.Pq())) return false;
  if (!Vq().isApprox(other.Vq())) return false;
  if (!Vv().isApprox(other.Vv())) return false;
  if (!Qdvdvff().isApprox(other.Qdvdvff())) return false;
  if (!Qfq().isApprox(other.Qfq())) return false;
  if (!Qqf().isApprox(other.Qqf())) return false;
  if (!Qxx().isApprox(other.Qxx())) return false;
  return true;
}


inline bool ImpulseKKTMatrix::hasNaN() const {
  if (Fqq_prev.hasNaN()) return true;
  if (F_.hasNaN()) return true;
  if (C_.hasNaN()) return true;
  if (Q_.hasNaN()) return true;
  return false;
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_KKT_MATRIX_HXX_ 