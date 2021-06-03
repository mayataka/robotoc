#ifndef IDOCP_IMPULSE_SPLIT_KKT_RESIDUAL_HXX_ 
#define IDOCP_IMPULSE_SPLIT_KKT_RESIDUAL_HXX_

#include "idocp/impulse/impulse_split_kkt_residual.hpp"


namespace idocp {

inline ImpulseSplitKKTResidual::ImpulseSplitKKTResidual(const Robot& robot) 
  : Fx(Eigen::VectorXd::Zero(2*robot.dimv())),
    lx(Eigen::VectorXd::Zero(2*robot.dimv())),
    ldv(Eigen::VectorXd::Zero(robot.dimv())),
    lf_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimi_(0) {
}


inline ImpulseSplitKKTResidual::ImpulseSplitKKTResidual() 
  : Fx(),
    lx(),
    ldv(),
    lf_full_(),
    dimv_(0), 
    dimi_(0) {
}


inline ImpulseSplitKKTResidual::~ImpulseSplitKKTResidual() {
}


inline void ImpulseSplitKKTResidual::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimi_  = impulse_status.dimf();
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::Fq() {
  assert(isDimensionConsistent());
  return Fx.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::Fq() const {
  assert(isDimensionConsistent());
  return Fx.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::Fv() {
  assert(isDimensionConsistent());
  return Fx.tail(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::Fv() const {
  assert(isDimensionConsistent());
  return Fx.tail(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::lq() {
  assert(isDimensionConsistent());
  return lx.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::lq() const {
  assert(isDimensionConsistent());
  return lx.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::lv() {
  assert(isDimensionConsistent());
  return lx.tail(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::lv() const {
  assert(isDimensionConsistent());
  return lx.tail(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitKKTResidual::lf() {
  return lf_full_.head(dimi_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitKKTResidual::lf() const {
  return lf_full_.head(dimi_);
}


inline void ImpulseSplitKKTResidual::setZero() {
  Fx.setZero();
  lx.setZero();
  ldv.setZero();
  lf().setZero();
}


inline int ImpulseSplitKKTResidual::dimi() const {
  return dimi_;
}


inline bool ImpulseSplitKKTResidual::isDimensionConsistent() const {
  if (Fx.size() != 2*dimv_) return false;
  if (lx.size() != 2*dimv_) return false;
  if (ldv.size() != dimv_) return false;
  return true;
}


inline bool ImpulseSplitKKTResidual::isApprox(
    const ImpulseSplitKKTResidual& other) const {
  assert(isDimensionConsistent());
  assert(other.isDimensionConsistent());
  if (!Fx.isApprox(other.Fx)) return false;
  if (!lx.isApprox(other.lx)) return false;
  if (!ldv.isApprox(other.ldv)) return false;
  if (!lf().isApprox(other.lf())) return false;
  return true;
}


inline bool ImpulseSplitKKTResidual::hasNaN() const {
  assert(isDimensionConsistent());
  if (Fx.hasNaN()) return true;
  if (lx.hasNaN()) return true;
  if (ldv.hasNaN()) return true;
  if (lf().hasNaN()) return true;
  return false;
}


inline void ImpulseSplitKKTResidual::setRandom() {
  Fx.setRandom();
  lx.setRandom();
  ldv.setRandom();
  lf().setRandom();
}


inline void ImpulseSplitKKTResidual::setRandom(
    const ImpulseStatus& impulse_status) {
  setImpulseStatus(impulse_status);
  setRandom();
}


inline ImpulseSplitKKTResidual ImpulseSplitKKTResidual::Random(
    const Robot& robot) {
  ImpulseSplitKKTResidual kkt_residual(robot);
  kkt_residual.setRandom();
  return kkt_residual;
}


inline ImpulseSplitKKTResidual ImpulseSplitKKTResidual::Random(
    const Robot& robot, const ImpulseStatus& impulse_status) {
  ImpulseSplitKKTResidual kkt_residual(robot);
  kkt_residual.setRandom(impulse_status);
  return kkt_residual;
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_SPLIT_KKT_RESIDUAL_HXX_ 