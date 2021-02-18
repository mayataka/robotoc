#ifndef IDOCP_SPLIT_DIRECTION_HXX_
#define IDOCP_SPLIT_DIRECTION_HXX_

#include "idocp/ocp/split_direction.hpp"

namespace idocp {

inline SplitDirection::SplitDirection(const Robot& robot) 
  : dnu_passive(Vector6d::Zero()),
    split_direction(Eigen::VectorXd::Zero(4*robot.dimv()+robot.dimu())),
    dxi_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    daf_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    dbetamu_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimu_(robot.dimu()), 
    dimx_(2*robot.dimv()), 
    dimf_(0), 
    dimi_(0),
    du_begin_(2*robot.dimv()), 
    dq_begin_(2*robot.dimv()+robot.dimu()), 
    dv_begin_(3*robot.dimv()+robot.dimu()),
    has_floating_base_(robot.hasFloatingBase()) {
}


inline SplitDirection::SplitDirection() 
  : dnu_passive(Vector6d::Zero()),
    split_direction(),
    dxi_full_(),
    daf_full_(),
    dbetamu_full_(),
    dimv_(0), 
    dimu_(0), 
    dimx_(0), 
    dimf_(0), 
    dimi_(0),
    du_begin_(0), 
    dq_begin_(0), 
    dv_begin_(0),
    has_floating_base_(false) {
}


inline SplitDirection::~SplitDirection() {
}


inline void SplitDirection::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
}


inline void SplitDirection::setContactStatusByDimension(const int dimf) {
  assert(dimf >= 0);
  assert(dimf % 3 == 0);
  dimf_ = dimf;
}


inline void SplitDirection::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimi_ = impulse_status.dimf();
}


inline void SplitDirection::setImpulseStatusByDimension(const int dimi) {
  assert(dimi >= 0);
  assert(dimi % 3 == 0);
  dimi_ = dimi;
}


inline void SplitDirection::setImpulseStatus() {
  dimi_ = 0;
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dlmdgmm() {
  return split_direction.head(dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dlmdgmm() const {
  return split_direction.head(dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dlmd() {
  return split_direction.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dlmd() const {
  return split_direction.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dgmm() {
  return split_direction.segment(dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dgmm() const {
  return split_direction.segment(dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::du() {
  return split_direction.segment(du_begin_, dimu_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::du() const {
  return split_direction.segment(du_begin_, dimu_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dq() {
  return split_direction.segment(dq_begin_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dq() const {
  return split_direction.segment(dq_begin_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dv() {
  return split_direction.segment(dv_begin_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dv() const {
  return split_direction.segment(dv_begin_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dx() {
  return split_direction.segment(dq_begin_, dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dx() const {
  return split_direction.segment(dq_begin_, dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dxi() {
  return dxi_full_.head(dimi_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dxi() const {
  return dxi_full_.head(dimi_);
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


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::df() {
  return daf_full_.segment(dimv_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::df() const {
  return daf_full_.segment(dimv_, dimf_);
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


inline void SplitDirection::setZero() {
  dnu_passive.setZero();
  split_direction.setZero();
  dxi_full_.setZero();
  daf_full_.setZero();
  dbetamu_full_.setZero();
}


inline int SplitDirection::dimf() const {
  return dimf_;
}


inline int SplitDirection::dimi() const {
  return dimi_;
}


inline bool SplitDirection::isApprox(const SplitDirection& other) const {
  assert(dimf()==other.dimf());
  assert(dimi()==other.dimi());
  if (!split_direction.isApprox(other.split_direction)) {
    return false;
  }
  if (!dxi().isApprox(other.dxi())) {
    return false;
  }
  if (!daf().isApprox(other.daf())) {
    return false;
  }
  if (!dbetamu().isApprox(other.dbetamu())) {
    return false;
  }
  if (has_floating_base_) {
    if (!dnu_passive.isApprox(other.dnu_passive)) {
      return false;
    }
  }
  return true;
}


inline void SplitDirection::setRandom() {
  split_direction.setRandom();
  dxi().setRandom();
  daf().setRandom();
  dbetamu().setRandom();
  if (has_floating_base_) {
    dnu_passive.setRandom();
  }
}


inline void SplitDirection::setRandom(const ContactStatus& contact_status) {
  setContactStatus(contact_status);
  setRandom();
}


inline void SplitDirection::setRandom(const ImpulseStatus& impulse_status) {
  setImpulseStatus(impulse_status);
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

} // namespace idocp 

#endif // IDOCP_SPLIT_OCP_DIRECTION_HXX_