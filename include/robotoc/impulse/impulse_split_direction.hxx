#ifndef ROBOTOC_IMPULSE_SPLIT_DIRECTION_HXX_
#define ROBOTOC_IMPULSE_SPLIT_DIRECTION_HXX_

#include "robotoc/impulse/impulse_split_direction.hpp"

#include <random>
#include <cassert>


namespace robotoc {

inline ImpulseSplitDirection::ImpulseSplitDirection(const Robot& robot) 
  : dx(Eigen::VectorXd::Zero(2*robot.dimv())),
    dlmdgmm(Eigen::VectorXd::Zero(2*robot.dimv())),
    dts(0.0),
    ddvf_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    dbetamu_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimi_(0) {
}


inline ImpulseSplitDirection::ImpulseSplitDirection() 
  : dx(),
    dlmdgmm(),
    dts(0.0),
    ddvf_full_(),
    dbetamu_full_(),
    dimv_(0), 
    dimi_(0) {
}


inline ImpulseSplitDirection::~ImpulseSplitDirection() {
}


inline void ImpulseSplitDirection::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimi_ = impulse_status.dimi();
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dq() {
  assert(isDimensionConsistent());
  return dx.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dq() const {
  assert(isDimensionConsistent());
  return dx.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dv() {
  assert(isDimensionConsistent());
  return dx.tail(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dv() const {
  assert(isDimensionConsistent());
  return dx.tail(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::ddvf() {
  assert(isDimensionConsistent());
  return ddvf_full_.head(dimv_+dimi_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::ddvf() const {
  return ddvf_full_.head(dimv_+dimi_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::ddv() {
  return ddvf_full_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::ddv() const {
  return ddvf_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::df() {
  return ddvf_full_.segment(dimv_, dimi_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::df() const {
  return ddvf_full_.segment(dimv_, dimi_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dlmd() {
  assert(isDimensionConsistent());
  return dlmdgmm.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dlmd() const {
  assert(isDimensionConsistent());
  return dlmdgmm.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dgmm() {
  assert(isDimensionConsistent());
  return dlmdgmm.tail(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dgmm() const {
  assert(isDimensionConsistent());
  return dlmdgmm.tail(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dbetamu() {
  return dbetamu_full_.head(dimv_+dimi_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dbetamu() const {
  return dbetamu_full_.head(dimv_+dimi_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dbeta() {
  return dbetamu_full_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dbeta() const {
  return dbetamu_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dmu() {
  return dbetamu_full_.segment(dimv_, dimi_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dmu() const {
  return dbetamu_full_.segment(dimv_, dimi_);
}


inline void ImpulseSplitDirection::setZero() {
  dx.setZero();
  ddvf().setZero();
  dlmdgmm.setZero();
  dbetamu().setZero();
}


inline int ImpulseSplitDirection::dimi() const {
  return dimi_;
}


inline bool ImpulseSplitDirection::isDimensionConsistent() const {
  if (dx.size() != 2*dimv_) return false;
  if (dlmdgmm.size() != 2*dimv_) return false;
  return true;
}


inline bool ImpulseSplitDirection::isApprox(
    const ImpulseSplitDirection& other) const {
  assert(isDimensionConsistent());
  assert(other.isDimensionConsistent());
  assert(dimi()==other.dimi());
  if (!dx.isApprox(other.dx)) return false;
  if (!ddvf().isApprox(other.ddvf())) return false;
  if (!dlmdgmm.isApprox(other.dlmdgmm)) return false;
  if (!dbetamu().isApprox(other.dbetamu())) return false;
  return true;
}


inline void ImpulseSplitDirection::setRandom() {
  assert(isDimensionConsistent());
  dx.setRandom();
  ddvf().setRandom();
  dlmdgmm.setRandom();
  dbetamu().setRandom();
}


inline void ImpulseSplitDirection::setRandom(
    const ImpulseStatus& impulse_status) {
  assert(isDimensionConsistent());
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

} // namespace robotoc 

#endif // ROBOTOC_IMPULSE_SPLIT_DIRECTION_HXX_ 