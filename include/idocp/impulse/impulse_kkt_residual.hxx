#ifndef IDOCP_IMPULSE_KKT_RESIDUAL_HXX_
#define IDOCP_IMPULSE_KKT_RESIDUAL_HXX_

#include "idocp/impulse/impulse_kkt_residual.hpp"

namespace idocp {

inline ImpulseKKTResidual::ImpulseKKTResidual(
    const Robot& robot, const bool use_contact_position_constraint) 
  : ldv(Eigen::VectorXd::Zero(robot.dimv())),
    dv_res(Eigen::VectorXd::Zero(robot.dimv())),
    kkt_residual_(Eigen::VectorXd::Zero(4*robot.dimv()+robot.max_dimf())),
    C_contact_velocity_full_(Eigen::VectorXd::Zero(robot.max_dimf())), 
    lf_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimf_(0), 
    dimc_(0),
    max_dimKKT_(4*robot.dimv()+robot.max_dimf()),
    dimKKT_(4*robot.dimv()),
    use_contact_position_constraint_(use_contact_position_constraint) {
}


inline ImpulseKKTResidual::ImpulseKKTResidual() 
  : ldv(),
    dv_res(),
    kkt_residual_(), 
    C_contact_velocity_full_(), 
    lf_full_(),
    dimv_(0), 
    dimx_(0), 
    dimf_(0), 
    dimc_(0),
    max_dimKKT_(0),
    dimKKT_(0),
    use_contact_position_constraint_(false) {
}


inline ImpulseKKTResidual::~ImpulseKKTResidual() {
}


inline void ImpulseKKTResidual::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
  if (use_contact_position_constraint_) {
    dimc_ = contact_status.dimf();
  }
  else {
    dimc_ = 0;
  }
  dimKKT_ = 4*dimv_ + dimc_;
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


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseKKTResidual::C_contact_position() {
  return kkt_residual_.segment(dimx_, dimc_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::lq() {
  return kkt_residual_.segment(dimx_+dimc_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::lv() {
  return kkt_residual_.segment(dimx_+dimc_+dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::lx() {
  return kkt_residual_.segment(dimx_+dimc_, dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseKKTResidual::C_contact_velocity() {
  return C_contact_velocity_full_.head(dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseKKTResidual::lf() {
  return lf_full_.head(dimf_);
}


inline void ImpulseKKTResidual::setZero() {
  ldv.setZero();
  dv_res.setZero();
  lf_full_.setZero();
  C_contact_velocity_full_.setZero();
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