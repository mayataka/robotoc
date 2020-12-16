#ifndef IDOCP_IMPULSE_SPLIT_KKT_RESIDUAL_HXX_ 
#define IDOCP_IMPULSE_SPLIT_KKT_RESIDUAL_HXX_

#include "idocp/impulse/impulse_split_kkt_residual.hpp"

namespace idocp {

inline ImpulseSplitKKTResidual::ImpulseSplitKKTResidual(const Robot& robot) 
  : KKT_residual(Eigen::VectorXd::Zero(4*robot.dimv())),
    ldv(Eigen::VectorXd::Zero(robot.dimv())),
    P_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    V_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    lf_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimf_(0), 
    dimKKT_(4*robot.dimv()) {
}


inline ImpulseSplitKKTResidual::ImpulseSplitKKTResidual() 
  : KKT_residual(),
    ldv(),
    P_full_(),
    V_full_(),
    lf_full_(),
    dimv_(0), 
    dimx_(0), 
    dimf_(0), 
    dimKKT_(0) {
}


inline ImpulseSplitKKTResidual::~ImpulseSplitKKTResidual() {
}


inline void ImpulseSplitKKTResidual::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimf_ = impulse_status.dimf();
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::Fq() {
  return KKT_residual.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::Fq() const {
  return KKT_residual.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::Fv() {
  return KKT_residual.segment(dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::Fv() const {
  return KKT_residual.segment(dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::Fx() {
  return KKT_residual.head(dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::Fx() const {
  return KKT_residual.head(dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::P() {
  return P_full_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::P() const {
  return P_full_.head(dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::V() {
  return V_full_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::V() const {
  return V_full_.head(dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::lq() {
  return KKT_residual.segment(dimx_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::lq() const {
  return KKT_residual.segment(dimx_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::lv() {
  return KKT_residual.segment(dimx_+dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::lv() const {
  return KKT_residual.segment(dimx_+dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::lx() {
  return KKT_residual.segment(dimx_, dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::lx() const {
  return KKT_residual.segment(dimx_, dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::lf() {
  return lf_full_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::lf() const {
  return lf_full_.head(dimf_);
}


inline void ImpulseSplitKKTResidual::setZero() {
  KKT_residual.setZero();
  ldv.setZero();
  lf_full_.setZero();
  P_full_.setZero();
  V_full_.setZero();
}


inline int ImpulseSplitKKTResidual::dimKKT() const {
  return dimKKT_;
}


inline int ImpulseSplitKKTResidual::dimf() const {
  return dimf_;
}


inline bool ImpulseSplitKKTResidual::isApprox(
    const ImpulseSplitKKTResidual& other) const {
  if (!Fx().isApprox(other.Fx())) return false;
  if (!lx().isApprox(other.lx())) return false;
  if (!ldv.isApprox(other.ldv)) return false;
  if (dimf_ > 0) {
    if (!P().isApprox(other.P())) return false;
    if (!V().isApprox(other.V())) return false;
    if (!lf().isApprox(other.lf())) return false;
  }
  return true;
}


inline bool ImpulseSplitKKTResidual::hasNaN() const {
  if (KKT_residual.hasNaN()) return true;
  if (ldv.hasNaN()) return true;
  if (P_full_.hasNaN()) return true;
  if (V_full_.hasNaN()) return true;
  if (lf_full_.hasNaN()) return true;
  return false;
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_SPLIT_KKT_RESIDUAL_HXX_ 