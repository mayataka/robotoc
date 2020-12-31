#ifndef IDOCP_SPLIT_UNKKT_RESIDUAL_HXX_ 
#define IDOCP_SPLIT_UNKKT_RESIDUAL_HXX_

#include "idocp/unocp/split_unkkt_residual.hpp"


namespace idocp {

inline SplitUnKKTResidual::SplitUnKKTResidual(const Robot& robot) 
  : KKT_residual(Eigen::VectorXd::Zero(5*robot.dimv())),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimKKT_(5*robot.dimv()) {
}


inline SplitUnKKTResidual::SplitUnKKTResidual() 
  : KKT_residual(),
    dimv_(0), 
    dimx_(0), 
    dimKKT_(0) {
}


inline SplitUnKKTResidual::~SplitUnKKTResidual() {
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitUnKKTResidual::Fq() {
  return KKT_residual.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitUnKKTResidual::Fq() const {
  return KKT_residual.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitUnKKTResidual::Fv() {
  return KKT_residual.segment(dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitUnKKTResidual::Fv() const {
  return KKT_residual.segment(dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitUnKKTResidual::Fx() {
  return KKT_residual.head(dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitUnKKTResidual::Fx() const {
  return KKT_residual.head(dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitUnKKTResidual::la() {
  return KKT_residual.segment(2*dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitUnKKTResidual::la() const {
  return KKT_residual.segment(2*dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitUnKKTResidual::lq() {
  return KKT_residual.segment(3*dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitUnKKTResidual::lq() const {
  return KKT_residual.segment(3*dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitUnKKTResidual::lv() {
  return KKT_residual.segment(4*dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitUnKKTResidual::lv() const {
  return KKT_residual.segment(4*dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitUnKKTResidual::lx() {
  return KKT_residual.segment(3*dimv_, dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitUnKKTResidual::lx() const {
  return KKT_residual.segment(3*dimv_, dimx_);
}


inline void SplitUnKKTResidual::setZero() {
  KKT_residual.setZero();
}


inline int SplitUnKKTResidual::dimKKT() const {
  return dimKKT_;
}


inline bool SplitUnKKTResidual::isApprox(
    const SplitUnKKTResidual& other) const {
  if (!KKT_residual.isApprox(other.KKT_residual)) return false;
  else return true;
}


inline bool SplitUnKKTResidual::hasNaN() const {
  if (KKT_residual.hasNaN()) return true;
  else return false;
}

} // namespace idocp 

#endif // IDOCP_SPLIT_UNKKT_RESIDUAL_HXX_ 