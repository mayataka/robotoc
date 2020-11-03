#ifndef IDOCP_IMPULSE_SPLIT_DIRECTION_HXX_
#define IDOCP_IMPULSE_SPLIT_DIRECTION_HXX_

#include "idocp/impulse/impulse_split_direction.hpp"

namespace idocp {

inline ImpulseSplitDirection::ImpulseSplitDirection(const Robot& robot) 
  : split_direction(Eigen::VectorXd::Zero(4*robot.dimv())),
    ddvf_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    dbetamu_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    dxi_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimf_(0), 
    dimKKT_(4*robot.dimv()),
    has_floating_base_(robot.has_floating_base()) {
}


inline ImpulseSplitDirection::ImpulseSplitDirection() 
  : split_direction(),
    ddvf_full_(),
    dbetamu_full_(),
    dxi_full_(),
    dimv_(0), 
    dimx_(0), 
    dimf_(0), 
    dimKKT_(0),
    has_floating_base_(false) {
}


inline ImpulseSplitDirection::~ImpulseSplitDirection() {
}


inline void ImpulseSplitDirection::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimf_ = impulse_status.dimp();
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dlmdgmm() {
  return split_direction.head(dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dlmdgmm() const {
  return split_direction.head(dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dlmd() {
  return split_direction.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dlmd() const {
  return split_direction.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dgmm() {
  return split_direction.segment(dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dgmm() const {
  return split_direction.segment(dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dq() {
  return split_direction.segment(dimx_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dq() const {
  return split_direction.segment(dimx_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dv() {
  return split_direction.segment(dimx_+dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dv() const {
  return split_direction.segment(dimx_+dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dx() {
  return split_direction.segment(dimx_, dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dx() const {
  return split_direction.segment(dimx_, dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::ddvf() {
  return ddvf_full_.head(dimv_+dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::ddvf() const {
  return ddvf_full_.head(dimv_+dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::ddv() {
  return ddvf_full_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::ddv() const {
  return ddvf_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::df() {
  return ddvf_full_.segment(dimv_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::df() const {
  return ddvf_full_.segment(dimv_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dbetamu() {
  return dbetamu_full_.head(dimv_+dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dbetamu() const {
  return dbetamu_full_.head(dimv_+dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dbeta() {
  return dbetamu_full_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dbeta() const {
  return dbetamu_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dmu() {
  return dbetamu_full_.segment(dimv_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dmu() const {
  return dbetamu_full_.segment(dimv_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dxi() {
  return dxi_full_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dxi() const {
  return dxi_full_.head(dimf_);
}


inline void ImpulseSplitDirection::setZero() {
  split_direction.setZero();
  ddvf_full_.setZero();
  dbetamu_full_.setZero();
  dxi_full_.setZero();
}


inline int ImpulseSplitDirection::dimKKT() const {
  return dimKKT_;
}


inline int ImpulseSplitDirection::dimf() const {
  return dimf_;
}


inline bool ImpulseSplitDirection::isApprox(
    const ImpulseSplitDirection& other) const {
  if (!dlmd().isApprox(other.dlmd())) {
    return false;
  }
  if (!dgmm().isApprox(other.dgmm())) {
    return false;
  }
  if (!dq().isApprox(other.dq())) {
    return false;
  }
  if (!dv().isApprox(other.dv())) {
    return false;
  }
  if (!ddv().isApprox(other.ddv())) {
    return false;
  }
  if (!df().isApprox(other.df())) {
    return false;
  }
  if (!dbeta().isApprox(other.dbeta())) {
    return false;
  }
  if (!dmu().isApprox(other.dmu())) {
    return false;
  }
  if (!dxi().isApprox(other.dxi())) {
    return false;
  }
  return true;
}


inline void ImpulseSplitDirection::setRandom() {
  split_direction.setRandom();
  ddvf().setRandom();
  dbetamu().setRandom();
  dxi().setRandom();
}


inline void ImpulseSplitDirection::setRandom(
    const ImpulseStatus& impulse_status) {
  setImpulseStatus(impulse_status);
  setRandom();
}


inline ImpulseSplitDirection ImpulseSplitDirection::Random(const Robot& robot) {
  ImpulseSplitDirection d(robot);
  d.setRandom();
  return d;
}


inline ImpulseSplitDirection ImpulseSplitDirection::Random(
    const Robot& robot, const ImpulseStatus& impulse_status) {
  ImpulseSplitDirection d(robot);
  d.setRandom(impulse_status);
  return d;
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_SPLIT_DIRECTION_HXX_ 