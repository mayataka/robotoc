#ifndef IDOCP_KKT_RESIDUAL_HXX_
#define IDOCP_KKT_RESIDUAL_HXX_

#include "idocp/ocp/kkt_residual.hpp"

namespace idocp {

inline KKTResidual::KKTResidual(const Robot& robot) 
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
    has_floating_base_(robot.has_floating_base()) {
}


inline KKTResidual::KKTResidual() 
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


inline KKTResidual::~KKTResidual() {
}


inline void KKTResidual::setContactStatus(const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
}


inline Eigen::VectorBlock<Eigen::VectorXd> KKTResidual::Fq() {
  return KKT_residual.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> KKTResidual::Fq() const {
  return KKT_residual.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> KKTResidual::Fv() {
  return KKT_residual.segment(dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> KKTResidual::Fv() const {
  return KKT_residual.segment(dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> KKTResidual::Fx() {
  return KKT_residual.head(dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> KKTResidual::Fx() const {
  return KKT_residual.head(dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> KKTResidual::lu() {
  return KKT_residual.segment(dimx_, dimu_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> KKTResidual::lu() const {
  return KKT_residual.segment(dimx_, dimu_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> KKTResidual::lq() {
  return KKT_residual.segment(dimx_+dimu_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> KKTResidual::lq() const {
  return KKT_residual.segment(dimx_+dimu_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> KKTResidual::lv() {
  return KKT_residual.segment(dimx_+dimu_+dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> KKTResidual::lv() const {
  return KKT_residual.segment(dimx_+dimu_+dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> KKTResidual::lx() {
  return KKT_residual.segment(dimx_+dimu_, dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> KKTResidual::lx() const {
  return KKT_residual.segment(dimx_+dimu_, dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> KKTResidual::lf() {
  return lf_full_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> KKTResidual::lf() const {
  return lf_full_.head(dimf_);
}


inline void KKTResidual::setZero() {
  KKT_residual.setZero();
  la.setZero();
  lu_passive.setZero();
  lf_full_.setZero();
}


inline int KKTResidual::dimKKT() const {
  return dimKKT_;
}


inline int KKTResidual::dimf() const {
  return dimf_;
}


inline bool KKTResidual::isApprox(const KKTResidual& other) const {
  if (!Fx().isApprox(other.Fx())) {
    return false;
  }
  if (!lu().isApprox(other.lu())) {
    return false;
  }
  if (!lx().isApprox(other.lx())) {
    return false;
  }
  if (!la.isApprox(other.la)) {
    return false;
  }
  if (!lf().isApprox(other.lf())) {
    return false;
  }
  if (has_floating_base_) {
    if (!lu_passive.isApprox(other.lu_passive)) {
      return false;
    }
  }
  return true;
}

} // namespace idocp 

#endif // IDOCP_KKT_RESIDUAL_HXX_