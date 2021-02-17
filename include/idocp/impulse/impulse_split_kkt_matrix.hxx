#ifndef IDOCP_IMPULSE_SPLIT_KKT_MATRIX_HXX_
#define IDOCP_IMPULSE_SPLIT_KKT_MATRIX_HXX_

#include "idocp/impulse/impulse_split_kkt_matrix.hpp"

#include <cassert>


namespace idocp {

inline ImpulseSplitKKTMatrix::ImpulseSplitKKTMatrix(const Robot& robot) 
  : Fqq_prev(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Fqq_inv(Matrix6d::Zero()),
    Fqq_prev_inv(Matrix6d::Zero()),
    FC_(Eigen::MatrixXd::Zero(2*robot.dimv()+robot.max_dimf(), 
                              2*robot.dimv()+robot.max_dimf())),
    Q_(Eigen::MatrixXd::Zero(3*robot.dimv()+robot.max_dimf(), 
                             3*robot.dimv()+robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimf_(0), 
    q_begin_(robot.dimv()),
    v_begin_(2*robot.dimv()),
    dimKKT_(4*robot.dimv()) {
}


inline ImpulseSplitKKTMatrix::ImpulseSplitKKTMatrix() 
  : Fqq_prev(),
    Fqq_inv(Matrix6d::Zero()),
    Fqq_prev_inv(Matrix6d::Zero()),
    FC_(),
    Q_(),
    dimv_(0), 
    dimx_(0), 
    dimf_(0), 
    q_begin_(0),
    v_begin_(0),
    dimKKT_(0) {
}


inline ImpulseSplitKKTMatrix::~ImpulseSplitKKTMatrix() {
}


inline void ImpulseSplitKKTMatrix::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimf_ = impulse_status.dimf();
  q_begin_ = dimv_ + dimf_;
  v_begin_ = 2*dimv_ + dimf_;
  dimKKT_ = 4*dimv_ + 2*dimf_;
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Fqf() {
  return FC_.block(0, 0, dimv_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Fqf() const {
  return FC_.block(0, 0, dimv_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Fqq() {
  return FC_.block(0, dimf_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Fqq() const {
  return FC_.block(0, dimf_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Fqv() {
  return FC_.block(0, dimf_+dimv_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Fqv() const {
  return FC_.block(0, dimf_+dimv_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Fvf() {
  return FC_.block(dimv_, 0, dimv_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Fvf() const {
  return FC_.block(dimv_, 0, dimv_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Fvq() {
  return FC_.block(dimv_, dimf_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Fvq() const {
  return FC_.block(dimv_, dimf_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Fvv() {
  return FC_.block(dimv_, dimf_+dimv_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Fvv() const {
  return FC_.block(dimv_, dimf_+dimv_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Fxf() {
  return FC_.block(0, 0, dimx_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Fxf() const {
  return FC_.block(0, 0, dimx_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Fxx() {
  return FC_.block(0, dimf_, dimx_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Fxx() const {
  return FC_.block(0, dimf_, dimx_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Vq() {
  return FC_.block(dimx_, dimf_, dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Vq() const {
  return FC_.block(dimx_, dimf_, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Vv() {
  return FC_.block(dimx_, dimf_+dimv_, dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Vv() const {
  return FC_.block(dimx_, dimf_+dimv_, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Qdvdvff() {
  return Q_.topLeftCorner(dimv_+dimf_, dimv_+dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Qdvdvff() const {
  return Q_.topLeftCorner(dimv_+dimf_, dimv_+dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Qdvdv() {
  return Q_.topLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Qdvdv() const {
  return Q_.topLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Qff() {
  return Q_.block(dimv_, dimv_, dimf_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Qff() const {
  return Q_.block(dimv_, dimv_, dimf_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Qfq() {
  return Q_.block(dimv_, q_begin_, dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Qfq() const {
  return Q_.block(dimv_, q_begin_, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Qfv() {
  return Q_.block(dimv_, v_begin_, dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Qfv() const {
  return Q_.block(dimv_, v_begin_, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Qqf() {
  return Q_.block(q_begin_, dimv_, dimv_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Qqf() const {
  return Q_.block(q_begin_, dimv_, dimv_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Qqq() {
  return Q_.block(q_begin_, q_begin_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Qqq() const {
  return Q_.block(q_begin_, q_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Qqv() {
  return Q_.block(q_begin_, v_begin_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Qqv() const {
  return Q_.block(q_begin_, v_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Qvf() {
  return Q_.block(v_begin_, dimv_, dimv_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Qvf() const {
  return Q_.block(v_begin_, dimv_, dimv_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Qvq() {
  return Q_.block(v_begin_, q_begin_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Qvq() const {
  return Q_.block(v_begin_, q_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Qvv() {
  return Q_.block(v_begin_, v_begin_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Qvv() const {
  return Q_.block(v_begin_, v_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Qxx() {
  return Q_.block(q_begin_, q_begin_, dimx_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Qxx() const {
  return Q_.block(q_begin_, q_begin_, dimx_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Qss() {
  return Q_.block(dimv_, dimv_, dimx_+dimf_, dimx_+dimf_); 
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Qss() const {
  return Q_.block(dimv_, dimv_, dimx_+dimf_, dimx_+dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseSplitKKTMatrix::Jac() {
  return FC_.topLeftCorner(dimx_+dimf_, dimx_+dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitKKTMatrix::Jac() const {
  return FC_.topLeftCorner(dimx_+dimf_, dimx_+dimf_); 
}


inline void ImpulseSplitKKTMatrix::symmetrize() {
  Q_.template triangularView<Eigen::StrictlyLower>() 
      = Q_.transpose().template triangularView<Eigen::StrictlyLower>();
}


inline void ImpulseSplitKKTMatrix::setZero() {
  Fqq_prev.setZero();
  FC_.setZero();
  Q_.setZero();
}


inline int ImpulseSplitKKTMatrix::dimKKT() const {
  return dimKKT_;
}


inline int ImpulseSplitKKTMatrix::dimf() const {
  return dimf_;
}


inline bool ImpulseSplitKKTMatrix::isApprox(
    const ImpulseSplitKKTMatrix& other) const {
  if (!Fxf().isApprox(other.Fxf())) return false;
  if (!Fxx().isApprox(other.Fxx())) return false;
  if (!Vq().isApprox(other.Vq())) return false;
  if (!Vv().isApprox(other.Vv())) return false;
  if (!Qdvdvff().isApprox(other.Qdvdvff())) return false;
  if (!Qfq().isApprox(other.Qfq())) return false;
  if (!Qqf().isApprox(other.Qqf())) return false;
  if (!Qxx().isApprox(other.Qxx())) return false;
  return true;
}


inline bool ImpulseSplitKKTMatrix::hasNaN() const {
  if (Fqq_prev.hasNaN()) return true;
  if (FC_.hasNaN()) return true;
  if (Q_.hasNaN()) return true;
  return false;
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_SPLIT_KKT_MATRIX_HXX_ 