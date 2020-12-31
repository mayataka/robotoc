#ifndef IDOCP_SPLIT_KKT_RESIDUAL_HXX_ 
#define IDOCP_SPLIT_KKT_RESIDUAL_HXX_

#include "idocp/ocp/split_kkt_residual.hpp"


namespace idocp {

inline SplitKKTResidual::SplitKKTResidual(const Robot& robot) 
  : KKT_residual(Eigen::VectorXd::Zero(4*robot.dimv()+robot.dimu())),
    la(Eigen::VectorXd::Zero(robot.dimv())),
    lu_passive(Vector6d::Zero()),
    lf_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimu_(robot.dimu()),
    dim_passive_(robot.dim_passive()),
    dimf_(0), 
    dimKKT_(4*robot.dimv()+robot.dimu()),
    has_floating_base_(robot.hasFloatingBase()) {
}


inline SplitKKTResidual::SplitKKTResidual() 
  : KKT_residual(),
    la(),
    lu_passive(Vector6d::Zero()),
    lf_full_(),
    dimv_(0), 
    dimx_(0), 
    dimu_(0),
    dim_passive_(0),
    dimf_(0), 
    dimKKT_(0),
    has_floating_base_(false) {
}


inline SplitKKTResidual::~SplitKKTResidual() {
}


inline void SplitKKTResidual::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::Fq() {
  return KKT_residual.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::Fq() const {
  return KKT_residual.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::Fv() {
  return KKT_residual.segment(dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::Fv() const {
  return KKT_residual.segment(dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::Fx() {
  return KKT_residual.head(dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::Fx() const {
  return KKT_residual.head(dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::lu() {
  return KKT_residual.segment(dimx_, dimu_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::lu() const {
  return KKT_residual.segment(dimx_, dimu_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::lq() {
  return KKT_residual.segment(dimx_+dimu_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::lq() const {
  return KKT_residual.segment(dimx_+dimu_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::lv() {
  return KKT_residual.segment(dimx_+dimu_+dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::lv() const {
  return KKT_residual.segment(dimx_+dimu_+dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::lx() {
  return KKT_residual.segment(dimx_+dimu_, dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::lx() const {
  return KKT_residual.segment(dimx_+dimu_, dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::lf() {
  return lf_full_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::lf() const {
  return lf_full_.head(dimf_);
}


inline void SplitKKTResidual::setZero() {
  KKT_residual.setZero();
  la.setZero();
  lu_passive.setZero();
  lf_full_.setZero();
}


inline int SplitKKTResidual::dimKKT() const {
  return dimKKT_;
}


inline int SplitKKTResidual::dimf() const {
  return dimf_;
}


inline bool SplitKKTResidual::isApprox(const SplitKKTResidual& other) const {
  if (!Fx().isApprox(other.Fx())) return false;
  if (!lu().isApprox(other.lu())) return false;
  if (!lx().isApprox(other.lx())) return false;
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
  if (KKT_residual.hasNaN()) return true;
  if (la.hasNaN()) return true;
  if (lu_passive.hasNaN()) return true;
  if (lf_full_.hasNaN()) return true;
  return false;
}

} // namespace idocp 

#endif // IDOCP_SPLIT_KKT_RESIDUAL_HXX_ 