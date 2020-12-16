#ifndef IDOCP_SPLIT_DIRECTION_HXX_
#define IDOCP_SPLIT_DIRECTION_HXX_

#include "idocp/ocp/split_direction.hpp"

namespace idocp {

inline SplitDirection::SplitDirection(const Robot& robot) 
  : split_direction(Eigen::VectorXd::Zero(4*robot.dimv()+robot.dimu())),
    du_passive(Vector6d::Zero()),
    dnu_passive(Vector6d::Zero()),
    daf_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    dbetamu_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimu_(robot.dimu()), 
    dimx_(2*robot.dimv()), 
    dimf_(0), 
    dimKKT_(4*robot.dimv()+robot.dimu()),
    has_floating_base_(robot.hasFloatingBase()) {
}


inline SplitDirection::SplitDirection() 
  : split_direction(),
    du_passive(Vector6d::Zero()),
    dnu_passive(Vector6d::Zero()),
    daf_full_(),
    dbetamu_full_(),
    dimv_(0), 
    dimu_(0), 
    dimx_(0), 
    dimf_(0), 
    dimKKT_(0),
    has_floating_base_(false) {
}


inline SplitDirection::~SplitDirection() {
}


inline void SplitDirection::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
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
  return split_direction.segment(dimx_, dimu_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::du() const {
  return split_direction.segment(dimx_, dimu_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dq() {
  return split_direction.segment(dimx_+dimu_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dq() const {
  return split_direction.segment(dimx_+dimu_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dv() {
  return split_direction.segment(dimx_+dimu_+dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dv() const {
  return split_direction.segment(dimx_+dimu_+dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dx() {
  return split_direction.segment(dimx_+dimu_, dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dx() const {
  return split_direction.segment(dimx_+dimu_, dimx_);
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
  split_direction.setZero();
  du_passive.setZero();
  dnu_passive.setZero();
  daf_full_.setZero();
  dbetamu_full_.setZero();
}


inline int SplitDirection::dimKKT() const {
  return dimKKT_;
}


inline int SplitDirection::dimf() const {
  return dimf_;
}


inline bool SplitDirection::isApprox(const SplitDirection& other) const {
  if (!dlmd().isApprox(other.dlmd())) {
    return false;
  }
  if (!dgmm().isApprox(other.dgmm())) {
    return false;
  }
  if (!du().isApprox(other.du())) {
    return false;
  }
  if (!dq().isApprox(other.dq())) {
    return false;
  }
  if (!dv().isApprox(other.dv())) {
    return false;
  }
  if (!da().isApprox(other.da())) {
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
  if (has_floating_base_) {
    if (!du_passive.isApprox(other.du_passive)) {
      return false;
    }
    if (!dnu_passive.isApprox(other.dnu_passive)) {
      return false;
    }
  }
  return true;
}


inline void SplitDirection::setRandom() {
  split_direction.setRandom();
  daf().setRandom();
  dbetamu().setRandom();
  if (has_floating_base_) {
    du_passive.setRandom();
    dnu_passive.setRandom();
  }
}


inline void SplitDirection::setRandom(const ContactStatus& contact_status) {
  setContactStatus(contact_status);
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

} // namespace idocp 

#endif // IDOCP_SPLIT_OCP_DIRECTION_HXX_