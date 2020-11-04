#ifndef IDOCP_IMPULSE_KKT_RESIDUAL_HXX_
#define IDOCP_IMPULSE_KKT_RESIDUAL_HXX_

#include "idocp/impulse/impulse_kkt_residual.hpp"

namespace idocp {

inline ImpulseKKTResidual::ImpulseKKTResidual(const Robot& robot) 
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


inline ImpulseKKTResidual::ImpulseKKTResidual() 
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


inline ImpulseKKTResidual::~ImpulseKKTResidual() {
}


inline void ImpulseKKTResidual::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimf_ = impulse_status.dimp();
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::Fq() {
  return KKT_residual.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseKKTResidual::Fq() const {
  return KKT_residual.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::Fv() {
  return KKT_residual.segment(dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseKKTResidual::Fv() const {
  return KKT_residual.segment(dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::Fx() {
  return KKT_residual.head(dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseKKTResidual::Fx() const {
  return KKT_residual.head(dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::P() {
  return P_full_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseKKTResidual::P() const {
  return P_full_.head(dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::V() {
  return V_full_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseKKTResidual::V() const {
  return V_full_.head(dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::lq() {
  return KKT_residual.segment(dimx_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseKKTResidual::lq() const {
  return KKT_residual.segment(dimx_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::lv() {
  return KKT_residual.segment(dimx_+dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseKKTResidual::lv() const {
  return KKT_residual.segment(dimx_+dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::lx() {
  return KKT_residual.segment(dimx_, dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseKKTResidual::lx() const {
  return KKT_residual.segment(dimx_, dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::lf() {
  return lf_full_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseKKTResidual::lf() const {
  return lf_full_.head(dimf_);
}


inline void ImpulseKKTResidual::setZero() {
  KKT_residual.setZero();
  ldv.setZero();
  lf_full_.setZero();
  P_full_.setZero();
  V_full_.setZero();
}


inline int ImpulseKKTResidual::dimKKT() const {
  return dimKKT_;
}


inline int ImpulseKKTResidual::dimf() const {
  return dimf_;
}


inline bool ImpulseKKTResidual::isApprox(
    const ImpulseKKTResidual& other) const {
  if (!Fx().isApprox(other.Fx())) {
    return false;
  }
  if (!lx().isApprox(other.lx())) {
    return false;
  }
  if (!ldv.isApprox(other.ldv)) {
    return false;
  }
  if (dimf_ > 0) {
    if (!P().isApprox(other.P())) {
      return false;
    }
    if (!V().isApprox(other.V())) {
      return false;
    }
    if (!lf().isApprox(other.lf())) {
      return false;
    }
  }
  return true;
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_KKT_RESIDUAL_HXX_ 