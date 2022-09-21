#ifndef ROBOTOC_SPLIT_DIRECTION_HXX_
#define ROBOTOC_SPLIT_DIRECTION_HXX_

#include "robotoc/ocp/split_direction.hpp"

#include <random>

namespace robotoc {

inline SplitDirection::SplitDirection(const Robot& robot) 
  : dx(Eigen::VectorXd::Zero(2*robot.dimv())),
    du(Eigen::VectorXd::Zero(robot.dimu())),
    dlmdgmm(Eigen::VectorXd::Zero(2*robot.dimv())),
    dnu_passive(Eigen::VectorXd::Zero(robot.dim_passive())),
    dts(0.0),
    dts_next(0.0),
    daf_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    dbetamu_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    dxi_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimu_(robot.dimu()), 
    dim_passive_(robot.dim_passive()), 
    dimf_(0), 
    dimi_(0) {
}


inline SplitDirection::SplitDirection() 
  : dx(),
    du(),
    dlmdgmm(),
    dnu_passive(),
    dts(0.0),
    dts_next(0.0),
    daf_full_(),
    dbetamu_full_(),
    dxi_full_(),
    dimv_(0), 
    dimu_(0), 
    dim_passive_(0), 
    dimf_(0), 
    dimi_(0) {
}


inline SplitDirection::~SplitDirection() {
}


inline void SplitDirection::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
}


inline void SplitDirection::setContactStatus(
    const ImpulseStatus& impulse_status) {
  dimf_ = impulse_status.dimi();
}


inline void SplitDirection::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimi_ = impulse_status.dimi();
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dq() {
  assert(isDimensionConsistent());
  return dx.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dq() const {
  assert(isDimensionConsistent());
  return dx.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dv() {
  assert(isDimensionConsistent());
  return dx.tail(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dv() const {
  assert(isDimensionConsistent());
  return dx.tail(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::daf() {
  return daf_full_.head(dimv_+dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::daf() const {
  return daf_full_.head(dimv_+dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::da() {
  return daf_full_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::da() const {
  return daf_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::ddvf() {
  return daf();
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::ddvf() const {
  return daf();
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::ddv() {
  return da();
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::ddv() const {
  return da();
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::df() {
  return daf_full_.segment(dimv_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::df() const {
  return daf_full_.segment(dimv_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dlmd() {
  assert(isDimensionConsistent());
  return dlmdgmm.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dlmd() const {
  assert(isDimensionConsistent());
  return dlmdgmm.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dgmm() {
  assert(isDimensionConsistent());
  return dlmdgmm.tail(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dgmm() const {
  assert(isDimensionConsistent());
  return dlmdgmm.tail(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dbetamu() {
  return dbetamu_full_.head(dimv_+dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dbetamu() const {
  return dbetamu_full_.head(dimv_+dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dbeta() {
  return dbetamu_full_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dbeta() const {
  return dbetamu_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dmu() {
  return dbetamu_full_.segment(dimv_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dmu() const {
  return dbetamu_full_.segment(dimv_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dxi() {
  return dxi_full_.head(dimi_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dxi() const {
  return dxi_full_.head(dimi_);
}


inline void SplitDirection::setZero() {
  dx.setZero();
  du.setZero();
  daf().setZero();
  dlmdgmm.setZero();
  dbetamu().setZero();
  dnu_passive.setZero();
  dxi().setZero();
  dts = 0.0;
  dts_next = 0.0;
}


inline int SplitDirection::dimf() const {
  return dimf_;
}


inline int SplitDirection::dimi() const {
  return dimi_;
}


inline bool SplitDirection::isDimensionConsistent() const {
  if (dx.size() != 2*dimv_) return false;
  if (du.size() != dimu_) return false;
  if (dlmdgmm.size() != 2*dimv_) return false;
  if (dnu_passive.size() != dim_passive_) return false;
  return true;
}


inline bool SplitDirection::isApprox(const SplitDirection& other) const {
  assert(isDimensionConsistent());
  assert(other.isDimensionConsistent());
  assert(dimf()==other.dimf());
  assert(dimi()==other.dimi());
  if (!dx.isApprox(other.dx)) {
    return false;
  }
  if (!du.isApprox(other.du)) {
    return false;
  }
  if (!daf().isApprox(other.daf())) {
    return false;
  }
  if (!dlmdgmm.isApprox(other.dlmdgmm)) {
    return false;
  }
  if (!dbetamu().isApprox(other.dbetamu())) {
    return false;
  }
  if (!dxi().isApprox(other.dxi())) {
    return false;
  }
  if (!dnu_passive.isApprox(other.dnu_passive)) {
    return false;
  }
  Eigen::VectorXd vec(2), other_vec(2);
  vec << dts, dts_next;
  other_vec << other.dts, other.dts_next;
  if (!vec.isApprox(other_vec)) {
    return false;
  } 
  return true;
}


inline void SplitDirection::setRandom() {
  assert(isDimensionConsistent());
  dx.setRandom();
  du.setRandom();
  daf().setRandom();
  dlmdgmm.setRandom();
  dbetamu().setRandom();
  dnu_passive.setRandom();
  dxi().setRandom();
  dts = Eigen::VectorXd::Random(1)[0];
  dts_next = Eigen::VectorXd::Random(1)[0];
}


inline void SplitDirection::setRandom(const ContactStatus& contact_status) {
  setContactStatus(contact_status);
  setRandom();
}


inline void SplitDirection::setRandom(const ImpulseStatus& impulse_status) {
  setContactStatus(impulse_status);
  setRandom();
}


inline void SplitDirection::setRandom(const ContactStatus& contact_status, 
                                      const ImpulseStatus& impulse_status) {
  setContactStatus(contact_status);
  setImpulseStatus(impulse_status);
  setRandom();
}


inline SplitDirection SplitDirection::Random(const Robot& robot) {
  SplitDirection d(robot);
  d.setRandom();
  return d;
}


inline SplitDirection SplitDirection::Random(
    const Robot& robot, const ContactStatus& contact_status) {
  SplitDirection d(robot);
  d.setRandom(contact_status);
  return d;
}


inline SplitDirection SplitDirection::Random(
    const Robot& robot, const ImpulseStatus& impulse_status) {
  SplitDirection d(robot);
  d.setRandom(impulse_status);
  return d;
}


inline SplitDirection SplitDirection::Random(
    const Robot& robot, const ContactStatus& contact_status, 
    const ImpulseStatus& impulse_status) {
  SplitDirection d(robot);
  d.setRandom(contact_status, impulse_status);
  return d;
}

} // namespace robotoc 

#endif // ROBOTOC_SPLIT_OCP_DIRECTION_HXX_