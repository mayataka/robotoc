#ifndef IDOCP_KKT_MATRIX_HXX_
#define IDOCP_KKT_MATRIX_HXX_

#include "idocp/ocp/kkt_matrix.hpp"

#include <cassert>


namespace idocp {

inline KKTMatrix::KKTMatrix(const Robot& robot) 
  : Fqq_prev(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    schur_complement_(2*robot.dimv(), 2*robot.dimv()+robot.dimu()),
    F_(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv()+robot.dimu())),
    C_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
    Q_(Eigen::MatrixXd::Zero(3*robot.dimv(), 3*robot.dimv())),
    Qaaff_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                      robot.dimv()+robot.max_dimf())),
    has_floating_base_(robot.has_floating_base()),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimu_(robot.dimu()), 
    dim_passive_(robot.dim_passive()),
    dimf_(0), 
    dimp_(0), 
    u_begin_(robot.dim_passive()),
    q_begin_(robot.dimv()),
    v_begin_(2*robot.dimv()),
    dimKKT_(4*robot.dimv()+robot.dimu()) {
}


inline KKTMatrix::KKTMatrix() 
  : Fqq_prev(),
    schur_complement_(),
    F_(),
    C_(),
    Q_(),
    Qaaff_full_(),
    has_floating_base_(false),
    dimv_(0), 
    dimx_(0), 
    dimu_(0), 
    dim_passive_(0),
    dimf_(0), 
    dimp_(0), 
    u_begin_(0),
    q_begin_(0),
    v_begin_(0),
    dimKKT_(0) {
}


inline KKTMatrix::~KKTMatrix() {
}


inline void KKTMatrix::setContactStatus(const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
}


inline void KKTMatrix::setImpulseStatus(const ImpulseStatus& impulse_status) {
  dimp_ = impulse_status.dimp();
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Fqu() {
  return F_.block(0, 0, dimv_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Fqu() const {
  return F_.block(0, 0, dimv_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Fqq() {
  return F_.block(0, dimu_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Fqq() const {
  return F_.block(0, dimu_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Fqv() {
  return F_.block(0, dimu_+dimv_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Fqv() const {
  return F_.block(0, dimu_+dimv_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Fvu() {
  return F_.block(dimv_, 0, dimv_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Fvu() const {
  return F_.block(dimv_, 0, dimv_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Fvq() {
  return F_.block(dimv_, dimu_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Fvq() const {
  return F_.block(dimv_, dimu_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Fvv() {
  return F_.block(dimv_, dimu_+dimv_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Fvv() const {
  return F_.block(dimv_, dimu_+dimv_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Fxu() {
  return F_.block(0, 0, dimx_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Fxu() const {
  return F_.block(0, 0, dimx_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Fxx() {
  return F_.block(0, dimu_, dimx_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Fxx() const {
  return F_.block(0, dimu_, dimx_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Cq() {
  return C_.topLeftCorner(dimp_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Cq() const {
  return C_.topLeftCorner(dimp_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Quu_full() {
  return Q_.topLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Quu_full() const {
  return Q_.topLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Quu_passive_topLeft() {
  return Q_.topLeftCorner(dim_passive_, dim_passive_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
KKTMatrix::Quu_passive_topLeft() const {
  return Q_.topLeftCorner(dim_passive_, dim_passive_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Quu_passive_topRight () {
  return Q_.block(0, dim_passive_, dim_passive_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
KKTMatrix::Quu_passive_topRight() const {
  return Q_.block(0, dim_passive_, dim_passive_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Quu_passive_bottomLeft() {
  return Q_.block(dim_passive_, 0, dimu_, dim_passive_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
KKTMatrix::Quu_passive_bottomLeft() const {
  return Q_.block(dim_passive_, 0, dimu_, dim_passive_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Quu() {
  return Q_.block(u_begin_, u_begin_, dimu_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Quu() const {
  return Q_.block(u_begin_, u_begin_, dimu_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Quq_full() {
  return Q_.block(0, q_begin_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Quq_full() const {
  return Q_.block(0, q_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Quq_passive() {
  return Q_.block(0, q_begin_, dim_passive_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
KKTMatrix::Quq_passive() const {
  return Q_.block(0, q_begin_, dim_passive_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Quq() {
  return Q_.block(u_begin_, q_begin_, dimu_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Quq() const {
  return Q_.block(u_begin_, q_begin_, dimu_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Quv_full() {
  return Q_.block(0, v_begin_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Quv_full() const {
  return Q_.block(0, v_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Quv_passive() {
  return Q_.block(0, v_begin_, dim_passive_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
KKTMatrix::Quv_passive() const {
  return Q_.block(0, v_begin_, dim_passive_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Quv() {
  return Q_.block(u_begin_, v_begin_, dimu_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Quv() const {
  return Q_.block(u_begin_, v_begin_, dimu_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Qqu_full() {
  return Q_.block(q_begin_, 0, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Qqu_full() const {
  return Q_.block(q_begin_, 0, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Qqu_passive() {
  return Q_.block(q_begin_, 0, dimv_, dim_passive_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
KKTMatrix::Qqu_passive() const {
  return Q_.block(q_begin_, 0, dimv_, dim_passive_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Qqu() {
  return Q_.block(q_begin_, u_begin_, dimv_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Qqu() const {
  return Q_.block(q_begin_, u_begin_, dimv_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Qqq() {
  return Q_.block(q_begin_, q_begin_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Qqq() const {
  return Q_.block(q_begin_, q_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Qqv() {
  return Q_.block(q_begin_, v_begin_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Qqv() const {
  return Q_.block(q_begin_, v_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Qvu_full() {
  return Q_.block(v_begin_, 0, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Qvu_full() const {
  return Q_.block(v_begin_, 0, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Qvu_passive() {
  return Q_.block(v_begin_, 0, dimv_, dim_passive_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
KKTMatrix::Qvu_passive() const {
  return Q_.block(v_begin_, 0, dimv_, dim_passive_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Qvu() {
  return Q_.block(v_begin_, u_begin_, dimv_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Qvu() const {
  return Q_.block(v_begin_, u_begin_, dimv_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Qvq() {
  return Q_.block(v_begin_, q_begin_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Qvq() const {
  return Q_.block(v_begin_, q_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Qvv() {
  return Q_.block(v_begin_, v_begin_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Qvv() const {
  return Q_.block(v_begin_, v_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Qux_full() {
  return Q_.block(0, q_begin_, dimv_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Qux_full() const {
  return Q_.block(0, q_begin_, dimv_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Qux_passive() {
  return Q_.block(0, q_begin_, dim_passive_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
KKTMatrix::Qux_passive() const {
  return Q_.block(0, q_begin_, dim_passive_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Qux() {
  return Q_.block(u_begin_, q_begin_, dimu_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Qux() const {
  return Q_.block(u_begin_, q_begin_, dimu_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Qxu_full() {
  return Q_.block(q_begin_, 0, dimx_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Qxu_full() const {
  return Q_.block(q_begin_, 0, dimx_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Qxu_passive() {
  return Q_.block(q_begin_, 0, dimx_, dim_passive_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
KKTMatrix::Qxu_passive() const {
  return Q_.block(q_begin_, 0, dimx_, dim_passive_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Qxu() {
  return Q_.block(q_begin_, u_begin_, dimx_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Qxu() const {
  return Q_.block(q_begin_, u_begin_, dimx_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Qxx() {
  return Q_.block(q_begin_, q_begin_, dimx_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Qxx() const {
  return Q_.block(q_begin_, q_begin_, dimx_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Qaaff() {
  return Qaaff_full_.topLeftCorner(dimv_+dimf_, dimv_+dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Qaaff() const {
  return Qaaff_full_.topLeftCorner(dimv_+dimf_, dimv_+dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Qaa() {
  return Qaaff_full_.topLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Qaa() const {
  return Qaaff_full_.topLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> KKTMatrix::Qff() {
  return Qaaff_full_.block(dimv_, dimv_, dimf_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> KKTMatrix::Qff() const {
  return Qaaff_full_.block(dimv_, dimv_, dimf_, dimf_);
}


inline void KKTMatrix::symmetrize() {
  Q_.template triangularView<Eigen::StrictlyLower>() 
      = Q_.transpose().template triangularView<Eigen::StrictlyLower>();
}


template <typename MatrixType>
inline void KKTMatrix::invert(
    const Eigen::MatrixBase<MatrixType>& KKT_matrix_inverse) {
  if (dimp_ > 0) {
  }
  else {
    assert(KKT_matrix_inverse.rows() == dimKKT_);
    assert(KKT_matrix_inverse.cols() == dimKKT_);
    schur_complement_.invertWithZeroTopLeftCorner(
        F_, Q_.bottomRightCorner(dimx_+dimu_, dimx_+dimu_), 
        const_cast<Eigen::MatrixBase<MatrixType>&>(KKT_matrix_inverse));
  }
}


inline void KKTMatrix::setZero() {
  Fqq_prev.setZero();
  F_.setZero();
  Q_.setZero();
  Qaaff_full_.setZero();
}


inline int KKTMatrix::dimKKT() const {
  return dimKKT_;
}


inline int KKTMatrix::dimf() const {
  return dimf_;
}


inline int KKTMatrix::dimp() const {
  return dimp_;
}


inline bool KKTMatrix::isApprox(const KKTMatrix& other) const {
  if (!Fxu().isApprox(other.Fxu())) {
    return false;
  }
  if (!Fxx().isApprox(other.Fxx())) {
    return false;
  }
  if (!Cq().isApprox(other.Cq())) {
    return false;
  }
  if (!Quu_full().isApprox(other.Quu_full())) {
    return false;
  }
  if (!Qux_full().isApprox(other.Qux_full())) {
    return false;
  }
  if (!Qxu_full().isApprox(other.Qxu_full())) {
    return false;
  }
  if (!Qxx().isApprox(other.Qxx())) {
    return false;
  }
  if (!Qaaff().isApprox(other.Qaaff())) {
    return false;
  }
  return true;
}

} // namespace idocp 

#endif // IDOCP_KKT_MATRIX_HXX_