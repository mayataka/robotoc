#ifndef IDOCP_SPLIT_KKT_MATRIX_HXX_
#define IDOCP_SPLIT_KKT_MATRIX_HXX_

#include "idocp/ocp/split_kkt_matrix.hpp"

#include <cassert>


namespace idocp {

inline SplitKKTMatrix::SplitKKTMatrix(const Robot& robot) 
  : Fqq_prev(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Fqq_prev_inv(
        Eigen::MatrixXd::Zero(robot.dim_passive(), robot.dim_passive())),
    F_(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv()+robot.dimu())),
    Pq_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
    Q_(Eigen::MatrixXd::Zero(3*robot.dimv(), 3*robot.dimv())),
    Qaaff_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                      robot.dimv()+robot.max_dimf())),
    has_floating_base_(robot.hasFloatingBase()),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimu_(robot.dimu()), 
    dim_passive_(robot.dim_passive()),
    dimf_(0), 
    dimi_(0), 
    dimKKT_(4*robot.dimv()+robot.dimu()),
    u_begin_(robot.dim_passive()),
    q_begin_(robot.dimv()),
    v_begin_(2*robot.dimv()) {
}


inline SplitKKTMatrix::SplitKKTMatrix() 
  : Fqq_prev(),
    Fqq_prev_inv(),
    F_(),
    Pq_full_(),
    Q_(),
    Qaaff_full_(),
    has_floating_base_(false),
    dimv_(0), 
    dimx_(0), 
    dimu_(0), 
    dim_passive_(0),
    dimf_(0), 
    dimi_(0), 
    dimKKT_(0),
    u_begin_(0),
    q_begin_(0),
    v_begin_(0) {
}


inline SplitKKTMatrix::~SplitKKTMatrix() {
}


inline void SplitKKTMatrix::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
}


inline void SplitKKTMatrix::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimi_ = impulse_status.dimf();
  dimKKT_ = 2*dimx_ + dimu_ + dimi_;
}


inline void SplitKKTMatrix::setImpulseStatus() {
  dimi_ = 0;
  dimKKT_ = 2*dimx_ + dimu_;
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Fqu() {
  return F_.block(0, 0, dimv_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Fqu() const {
  return F_.block(0, 0, dimv_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Fqq() {
  return F_.block(0, dimu_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Fqq() const {
  return F_.block(0, dimu_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Fqv() {
  return F_.block(0, dimu_+dimv_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Fqv() const {
  return F_.block(0, dimu_+dimv_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Fvu() {
  return F_.block(dimv_, 0, dimv_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Fvu() const {
  return F_.block(dimv_, 0, dimv_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Fvq() {
  return F_.block(dimv_, dimu_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Fvq() const {
  return F_.block(dimv_, dimu_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Fvv() {
  return F_.block(dimv_, dimu_+dimv_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Fvv() const {
  return F_.block(dimv_, dimu_+dimv_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Fxu() {
  return F_.block(0, 0, dimx_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Fxu() const {
  return F_.block(0, 0, dimx_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Fxx() {
  return F_.block(0, dimu_, dimx_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Fxx() const {
  return F_.block(0, dimu_, dimx_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Pq() {
  return Pq_full_.topLeftCorner(dimi_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Pq() const {
  return Pq_full_.topLeftCorner(dimi_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Quu_full() {
  return Q_.topLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitKKTMatrix::Quu_full() const {
  return Q_.topLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Quu_passive_topLeft() {
  return Q_.topLeftCorner(dim_passive_, dim_passive_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitKKTMatrix::Quu_passive_topLeft() const {
  return Q_.topLeftCorner(dim_passive_, dim_passive_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Quu_passive_topRight () {
  return Q_.block(0, dim_passive_, dim_passive_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitKKTMatrix::Quu_passive_topRight() const {
  return Q_.block(0, dim_passive_, dim_passive_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Quu_passive_bottomLeft() {
  return Q_.block(dim_passive_, 0, dimu_, dim_passive_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitKKTMatrix::Quu_passive_bottomLeft() const {
  return Q_.block(dim_passive_, 0, dimu_, dim_passive_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Quu() {
  return Q_.block(u_begin_, u_begin_, dimu_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Quu() const {
  return Q_.block(u_begin_, u_begin_, dimu_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Quq_full() {
  return Q_.block(0, q_begin_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitKKTMatrix::Quq_full() const {
  return Q_.block(0, q_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Quq_passive() {
  return Q_.block(0, q_begin_, dim_passive_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitKKTMatrix::Quq_passive() const {
  return Q_.block(0, q_begin_, dim_passive_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Quq() {
  return Q_.block(u_begin_, q_begin_, dimu_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Quq() const {
  return Q_.block(u_begin_, q_begin_, dimu_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Quv_full() {
  return Q_.block(0, v_begin_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitKKTMatrix::Quv_full() const {
  return Q_.block(0, v_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Quv_passive() {
  return Q_.block(0, v_begin_, dim_passive_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitKKTMatrix::Quv_passive() const {
  return Q_.block(0, v_begin_, dim_passive_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Quv() {
  return Q_.block(u_begin_, v_begin_, dimu_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Quv() const {
  return Q_.block(u_begin_, v_begin_, dimu_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qqu_full() {
  return Q_.block(q_begin_, 0, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitKKTMatrix::Qqu_full() const {
  return Q_.block(q_begin_, 0, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qqu_passive() {
  return Q_.block(q_begin_, 0, dimv_, dim_passive_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitKKTMatrix::Qqu_passive() const {
  return Q_.block(q_begin_, 0, dimv_, dim_passive_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qqu() {
  return Q_.block(q_begin_, u_begin_, dimv_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qqu() const {
  return Q_.block(q_begin_, u_begin_, dimv_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qqq() {
  return Q_.block(q_begin_, q_begin_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qqq() const {
  return Q_.block(q_begin_, q_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qqv() {
  return Q_.block(q_begin_, v_begin_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qqv() const {
  return Q_.block(q_begin_, v_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qvu_full() {
  return Q_.block(v_begin_, 0, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitKKTMatrix::Qvu_full() const {
  return Q_.block(v_begin_, 0, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qvu_passive() {
  return Q_.block(v_begin_, 0, dimv_, dim_passive_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitKKTMatrix::Qvu_passive() const {
  return Q_.block(v_begin_, 0, dimv_, dim_passive_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qvu() {
  return Q_.block(v_begin_, u_begin_, dimv_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qvu() const {
  return Q_.block(v_begin_, u_begin_, dimv_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qvq() {
  return Q_.block(v_begin_, q_begin_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qvq() const {
  return Q_.block(v_begin_, q_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qvv() {
  return Q_.block(v_begin_, v_begin_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qvv() const {
  return Q_.block(v_begin_, v_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qux_full() {
  return Q_.block(0, q_begin_, dimv_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitKKTMatrix::Qux_full() const {
  return Q_.block(0, q_begin_, dimv_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qux_passive() {
  return Q_.block(0, q_begin_, dim_passive_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitKKTMatrix::Qux_passive() const {
  return Q_.block(0, q_begin_, dim_passive_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qux() {
  return Q_.block(u_begin_, q_begin_, dimu_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qux() const {
  return Q_.block(u_begin_, q_begin_, dimu_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qxu_full() {
  return Q_.block(q_begin_, 0, dimx_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitKKTMatrix::Qxu_full() const {
  return Q_.block(q_begin_, 0, dimx_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qxu_passive() {
  return Q_.block(q_begin_, 0, dimx_, dim_passive_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitKKTMatrix::Qxu_passive() const {
  return Q_.block(q_begin_, 0, dimx_, dim_passive_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qxu() {
  return Q_.block(q_begin_, u_begin_, dimx_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qxu() const {
  return Q_.block(q_begin_, u_begin_, dimx_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qxx() {
  return Q_.block(q_begin_, q_begin_, dimx_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qxx() const {
  return Q_.block(q_begin_, q_begin_, dimx_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qaaff() {
  return Qaaff_full_.topLeftCorner(dimv_+dimf_, dimv_+dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qaaff() const {
  return Qaaff_full_.topLeftCorner(dimv_+dimf_, dimv_+dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qaa() {
  return Qaaff_full_.topLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qaa() const {
  return Qaaff_full_.topLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qff() {
  return Qaaff_full_.block(dimv_, dimv_, dimf_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qff() const {
  return Qaaff_full_.block(dimv_, dimv_, dimf_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qss() {
  return Q_.bottomRightCorner(dimx_+dimu_, dimx_+dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qss() const {
  return Q_.bottomRightCorner(dimx_+dimu_, dimx_+dimu_);
}


inline Eigen::MatrixXd& SplitKKTMatrix::Jac() {
  return F_;
}


inline const Eigen::MatrixXd& SplitKKTMatrix::Jac() const {
  return F_;
}


inline void SplitKKTMatrix::symmetrize() {
  Q_.template triangularView<Eigen::StrictlyLower>() 
      = Q_.transpose().template triangularView<Eigen::StrictlyLower>();
}


inline void SplitKKTMatrix::setZero() {
  Fqq_prev.setZero();
  F_.setZero();
  Pq_full_.setZero();
  Q_.setZero();
  Qaaff_full_.setZero();
}


inline int SplitKKTMatrix::dimKKT() const {
  return dimKKT_;
}


inline int SplitKKTMatrix::dimf() const {
  return dimf_;
}


inline int SplitKKTMatrix::dimi() const {
  return dimi_;
}


inline bool SplitKKTMatrix::isApprox(const SplitKKTMatrix& other) const {
  if (!Fqq_prev.isApprox(other.Fqq_prev)) return false;
  if (!Jac().isApprox(other.Jac())) return false;
  if (!Pq().isApprox(other.Pq())) return false;
  if (!Quu_full().isApprox(other.Quu_full())) return false;
  if (!Qux_full().isApprox(other.Qux_full())) return false;
  if (!Qxu_full().isApprox(other.Qxu_full())) return false;
  if (!Qxx().isApprox(other.Qxx())) return false;
  if (!Qaaff().isApprox(other.Qaaff())) return false;
  return true;
}


inline bool SplitKKTMatrix::hasNaN() const {
  if (Fqq_prev.hasNaN()) return true;
  if (F_.hasNaN()) return true;
  if (Pq_full_.hasNaN()) return true;
  if (Q_.hasNaN()) return true;
  if (Qaaff_full_.hasNaN()) return true;
  return false;
}

} // namespace idocp 

#endif // IDOCP_SPLIT_KKT_MATRIX_HXX_ 