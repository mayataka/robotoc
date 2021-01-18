#ifndef IDOCP_IMPULSE_SPLIT_KKT_RESIDUAL_HXX_ 
#define IDOCP_IMPULSE_SPLIT_KKT_RESIDUAL_HXX_

#include "idocp/impulse/impulse_split_kkt_residual.hpp"


namespace idocp {

inline ImpulseSplitKKTResidual::ImpulseSplitKKTResidual(const Robot& robot) 
  : ldv(Eigen::VectorXd::Zero(robot.dimv())),
    KKT_residual_full_(Eigen::VectorXd::Zero(4*robot.dimv()+2*robot.max_dimf())),
    P_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimf_(0), 
    dimKKT_(4*robot.dimv()), 
    lf_begin_(dimx_), 
    lq_begin_(dimx_), 
    lv_begin_(dimx_+dimv_) {
}


inline ImpulseSplitKKTResidual::ImpulseSplitKKTResidual() 
  : ldv(),
    KKT_residual_full_(),
    P_full_(),
    dimv_(0), 
    dimx_(0), 
    dimf_(0), 
    dimKKT_(0),
    lf_begin_(0), 
    lq_begin_(0), 
    lv_begin_(0) {
}


inline ImpulseSplitKKTResidual::~ImpulseSplitKKTResidual() {
}


inline void ImpulseSplitKKTResidual::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimf_     = impulse_status.dimf();
  dimKKT_   = 2*dimx_ + 2*dimf_;
  lf_begin_ = dimx_ + dimf_;
  lq_begin_ = dimx_ + 2*dimf_;
  lv_begin_ = dimx_ + 2*dimf_ + dimv_;
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::Fq() {
  return KKT_residual_full_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::Fq() const {
  return KKT_residual_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::Fv() {
  return KKT_residual_full_.segment(dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::Fv() const {
  return KKT_residual_full_.segment(dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::Fx() {
  return KKT_residual_full_.head(dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::Fx() const {
  return KKT_residual_full_.head(dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::V() {
  return KKT_residual_full_.segment(dimx_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::V() const {
  return KKT_residual_full_.segment(dimx_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::lf() {
  return KKT_residual_full_.segment(lf_begin_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::lf() const {
  return KKT_residual_full_.segment(lf_begin_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::lq() {
  return KKT_residual_full_.segment(lq_begin_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::lq() const {
  return KKT_residual_full_.segment(lq_begin_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::lv() {
  return KKT_residual_full_.segment(lv_begin_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::lv() const {
  return KKT_residual_full_.segment(lv_begin_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::lx() {
  return KKT_residual_full_.segment(lq_begin_, dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::lx() const {
  return KKT_residual_full_.segment(lq_begin_, dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::P() {
  return P_full_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::P() const {
  return P_full_.head(dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseSplitKKTResidual::KKT_residual() {
  return KKT_residual_full_.head(dimKKT_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::KKT_residual() const {
  return KKT_residual_full_.head(dimKKT_);
}


inline void ImpulseSplitKKTResidual::setZero() {
  KKT_residual_full_.setZero();
  ldv.setZero();
  P_full_.setZero();
}


inline int ImpulseSplitKKTResidual::dimKKT() const {
  return dimKKT_;
}


inline int ImpulseSplitKKTResidual::dimf() const {
  return dimf_;
}


inline bool ImpulseSplitKKTResidual::isApprox(
    const ImpulseSplitKKTResidual& other) const {
  if (!KKT_residual().isApprox(other.KKT_residual())) return false;
  if (!ldv.isApprox(other.ldv)) return false;
  if (dimf_ > 0) {
    if (!P().isApprox(other.P())) return false;
  }
  return true;
}


inline bool ImpulseSplitKKTResidual::hasNaN() const {
  if (KKT_residual_full_.hasNaN()) return true;
  if (ldv.hasNaN()) return true;
  if (P_full_.hasNaN()) return true;
  return false;
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_SPLIT_KKT_RESIDUAL_HXX_ 