#ifndef IDOCP_KKT_RESIDUAL_HXX_
#define IDOCP_KKT_RESIDUAL_HXX_

#include "idocp/ocp/kkt_residual.hpp"

namespace idocp {

inline KKTResidual::KKTResidual(const Robot& robot) 
  : lu(Eigen::VectorXd::Zero(robot.dimv())),
    u_res(Eigen::VectorXd::Zero(robot.dimv())),
    kkt_residual_(Eigen::VectorXd::Zero(
        5*robot.dimv()+robot.dim_passive()+2*robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimf_(robot.dimf()), 
    dimc_(robot.dim_passive()+robot.dimf()),
    max_dimKKT_(5*robot.dimv()+robot.dim_passive()+2*robot.max_dimf()),
    dimKKT_(5*robot.dimv()+robot.dim_passive()+2*robot.dimf()) {
}


inline KKTResidual::KKTResidual() 
  : lu(),
    u_res(),
    kkt_residual_(), 
    dimv_(0), 
    dimx_(0), 
    dimf_(0), 
    dimc_(0),
    max_dimKKT_(0),
    dimKKT_(0) {
}


inline KKTResidual::~KKTResidual() {
}


inline void KKTResidual::setContactStatus(const Robot& robot) {
  dimf_ = robot.dimf();
  dimc_ = robot.dim_passive() + robot.dimf();
  dimKKT_ = 5*robot.dimv() + robot.dim_passive() + 2*robot.dimf();
}


inline Eigen::VectorBlock<Eigen::VectorXd> KKTResidual::KKT_residual() {
  return kkt_residual_.head(dimKKT_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> KKTResidual::Fq() {
  return kkt_residual_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> KKTResidual::Fv() {
  return kkt_residual_.segment(dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> KKTResidual::Fx() {
  return kkt_residual_.segment(0, dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> KKTResidual::C() {
  return kkt_residual_.segment(dimx_, dimc_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> KKTResidual::la() {
  return kkt_residual_.segment(dimx_+dimc_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> KKTResidual::lf() {
  return kkt_residual_.segment(dimx_+dimc_+dimv_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> KKTResidual::lq() {
  return kkt_residual_.segment(dimx_+dimc_+dimv_+dimf_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> KKTResidual::lv() {
  return kkt_residual_.segment(dimx_+dimc_+2*dimv_+dimf_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> KKTResidual::lx() {
  return kkt_residual_.segment(dimx_+dimc_+dimv_+dimf_, dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> KKTResidual::laf() {
  return kkt_residual_.segment(dimx_+dimc_, dimv_+dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
KKTResidual::KKT_residual() const {
  return kkt_residual_.head(dimKKT_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
KKTResidual::Fq() const {
  return kkt_residual_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
KKTResidual::Fv() const {
  return kkt_residual_.segment(dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
KKTResidual::Fx() const {
  return kkt_residual_.segment(0, dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
KKTResidual::C() const {
  return kkt_residual_.segment(dimx_, dimc_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
KKTResidual::la() const {
  return kkt_residual_.segment(dimx_+dimc_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
KKTResidual::lf() const {
  return kkt_residual_.segment(dimx_+dimc_+dimv_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
KKTResidual::lq() const {
  return kkt_residual_.segment(dimx_+dimc_+dimv_+dimf_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
KKTResidual::lv() const {
  return kkt_residual_.segment(dimx_+dimc_+2*dimv_+dimf_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
KKTResidual::lx() const {
  return kkt_residual_.segment(dimx_+dimc_+dimv_+dimf_, dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
KKTResidual::laf() const {
  return kkt_residual_.segment(dimx_+dimc_, dimv_+dimf_);
}


inline double KKTResidual::squaredKKTErrorNorm(const double dtau) const {
  assert(dtau > 0);
  double error = kkt_residual_.head(dimKKT_).squaredNorm();
  error += lu.squaredNorm();
  error += dtau * dtau * u_res.squaredNorm();
  return error;
}


inline void KKTResidual::setZeroMinimum() {
  lu.setZero();
  kkt_residual_.segment(dimx_+dimc_, 3*dimv_+dimf_).setZero();
}


inline void KKTResidual::setZero() {
  lu.setZero();
  u_res.setZero();
  kkt_residual_.setZero();
}


inline int KKTResidual::dimKKT() const {
  return dimKKT_;
}


inline int KKTResidual::max_dimKKT() const {
  return max_dimKKT_;
}


inline int KKTResidual::dimc() const {
  return dimc_;
}


inline int KKTResidual::dimf() const {
  return dimf_;
}

} // namespace idocp 

#endif // IDOCP_KKT_RESIDUAL_HXX_