#ifndef IDOCP_SPLIT_KKT_RESIDUAL_HXX_ 
#define IDOCP_SPLIT_KKT_RESIDUAL_HXX_

#include "idocp/ocp/split_kkt_residual.hpp"


namespace idocp {

inline SplitKKTResidual::SplitKKTResidual(const Robot& robot) 
  : la(Eigen::VectorXd::Zero(robot.dimv())),
    lu_passive(Vector6d::Zero()),
    kkt_residual_full_(
        Eigen::VectorXd::Zero(4*robot.dimv()+robot.dimu()+robot.max_dimf())),
    lf_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimu_(robot.dimu()),
    dim_passive_(robot.dim_passive()),
    dimf_(0), 
    dimi_(0), 
    dimKKT_(4*robot.dimv()+robot.dimu()),
    lu_begin_(2*robot.dimv()), 
    lq_begin_(2*robot.dimv()+robot.dimu()), 
    lv_begin_(3*robot.dimv()+robot.dimu()),
    has_floating_base_(robot.hasFloatingBase()) {
}


inline SplitKKTResidual::SplitKKTResidual() 
  : la(),
    lu_passive(Vector6d::Zero()),
    lf_full_(),
    dimv_(0), 
    dimx_(0), 
    dimu_(0),
    dim_passive_(0),
    dimf_(0), 
    dimi_(0), 
    dimKKT_(0),
    lu_begin_(0), 
    lq_begin_(0), 
    lv_begin_(0),
    has_floating_base_(false) {
}


inline SplitKKTResidual::~SplitKKTResidual() {
}


inline void SplitKKTResidual::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
}


inline void SplitKKTResidual::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimi_ = impulse_status.dimf();
  dimKKT_ = 2*dimx_ + dimu_ + dimi_;
  lu_begin_ = dimx_ + dimi_;
  lq_begin_ = dimx_ + dimi_ + dimu_; 
  lv_begin_ = dimx_ + dimi_ + dimu_ + dimv_; 
}


inline void SplitKKTResidual::setImpulseStatus() {
  dimi_ = 0;
  dimKKT_ = 2*dimx_ + dimu_;
  lu_begin_ = dimx_;
  lq_begin_ = dimx_ + dimu_; 
  lv_begin_ = dimx_ + dimu_ + dimv_; 
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::Fq() {
  return kkt_residual_full_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::Fq() const {
  return kkt_residual_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::Fv() {
  return kkt_residual_full_.segment(dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::Fv() const {
  return kkt_residual_full_.segment(dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::Fx() {
  return kkt_residual_full_.head(dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::Fx() const {
  return kkt_residual_full_.head(dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::P() {
  return kkt_residual_full_.segment(dimx_, dimi_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::P() const {
  return kkt_residual_full_.segment(dimx_, dimi_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::lu() {
  return kkt_residual_full_.segment(lu_begin_, dimu_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::lu() const {
  return kkt_residual_full_.segment(lu_begin_, dimu_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::lq() {
  return kkt_residual_full_.segment(lq_begin_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::lq() const {
  return kkt_residual_full_.segment(lq_begin_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::lv() {
  return kkt_residual_full_.segment(lv_begin_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::lv() const {
  return kkt_residual_full_.segment(lv_begin_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::lx() {
  return kkt_residual_full_.segment(lq_begin_, dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::lx() const {
  return kkt_residual_full_.segment(lq_begin_, dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
SplitKKTResidual::splitKKTResidual() {
  return kkt_residual_full_.head(dimKKT_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::splitKKTResidual() const {
  return kkt_residual_full_.head(dimKKT_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::lf() {
  return lf_full_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::lf() const {
  return lf_full_.head(dimf_);
}


inline void SplitKKTResidual::setZero() {
  la.setZero();
  lu_passive.setZero();
  kkt_residual_full_.setZero();
  lf_full_.setZero();
}


inline int SplitKKTResidual::dimKKT() const {
  return dimKKT_;
}


inline int SplitKKTResidual::dimf() const {
  return dimf_;
}


inline int SplitKKTResidual::dimi() const {
  return dimi_;
}


inline bool SplitKKTResidual::isApprox(const SplitKKTResidual& other) const {
  if (!splitKKTResidual().isApprox(other.splitKKTResidual())) return false;
  if (!la.isApprox(other.la)) return false;
  if (dimf_ > 0) {
    if (!lf().isApprox(other.lf())) return false;
  }
  if (has_floating_base_) {
    if (!lu_passive.isApprox(other.lu_passive)) return false;
  }
  return true;
}


inline bool SplitKKTResidual::hasNaN() const {
  if (la.hasNaN()) return true;
  if (lu_passive.hasNaN()) return true;
  if (kkt_residual_full_.hasNaN()) return true;
  if (lf_full_.hasNaN()) return true;
  return false;
}

} // namespace idocp 

#endif // IDOCP_SPLIT_KKT_RESIDUAL_HXX_ 