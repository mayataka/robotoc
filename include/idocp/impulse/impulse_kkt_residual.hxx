#ifndef IDOCP_IMPULSE_KKT_RESIDUAL_HXX_
#define IDOCP_IMPULSE_KKT_RESIDUAL_HXX_

#include "idocp/impulse/impulse_kkt_residual.hpp"

namespace idocp {

inline ImpulseKKTResidual::ImpulseKKTResidual(const Robot& robot) 
  : ldv(Eigen::VectorXd::Zero(robot.dimv())),
    dv_res(Eigen::VectorXd::Zero(robot.dimv())),
    kkt_residual_(Eigen::VectorXd::Zero(4*robot.dimv()+2*robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimf_(0), 
    dimc_(0),
    max_dimKKT_(4*robot.dimv()+2*robot.max_dimf()),
    dimKKT_(4*robot.dimv()) {
}


inline ImpulseKKTResidual::ImpulseKKTResidual() 
  : ldv(),
    dv_res(),
    kkt_residual_(), 
    dimv_(0), 
    dimx_(0), 
    dimf_(0), 
    dimc_(0),
    max_dimKKT_(0),
    dimKKT_(0) {
}


inline ImpulseKKTResidual::~ImpulseKKTResidual() {
}


inline void ImpulseKKTResidual::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
  dimc_ = contact_status.dimf();
  dimKKT_ = 4*dimv_ + dimf_ + dimc_;
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::KKT_residual() {
  return kkt_residual_.head(dimKKT_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::Fq() {
  return kkt_residual_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::Fv() {
  return kkt_residual_.segment(dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::Fx() {
  return kkt_residual_.segment(0, dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseKKTResidual::Fx() const {
  return kkt_residual_.segment(0, dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::C() {
  return kkt_residual_.segment(dimx_, dimc_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseKKTResidual::C() const {
  return kkt_residual_.segment(dimx_, dimc_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::lf() {
  return kkt_residual_.segment(dimx_+dimc_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::lq() {
  return kkt_residual_.segment(dimx_+dimc_+dimf_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::lv() {
  return kkt_residual_.segment(dimx_+dimc_+dimf_+dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::lx() {
  return kkt_residual_.segment(dimx_+dimc_+dimf_, dimx_);
}


inline void ImpulseKKTResidual::setZero() {
  ldv.setZero();
  dv_res.setZero();
  kkt_residual_.setZero();
}


inline int ImpulseKKTResidual::dimKKT() const {
  return dimKKT_;
}


inline int ImpulseKKTResidual::max_dimKKT() const {
  return max_dimKKT_;
}


inline int ImpulseKKTResidual::dimc() const {
  return dimc_;
}


inline int ImpulseKKTResidual::dimf() const {
  return dimf_;
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_KKT_RESIDUAL_HXX_ 