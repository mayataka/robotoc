#ifndef IDOCP_SPLIT_DIRECTION_HXX_
#define IDOCP_SPLIT_DIRECTION_HXX_

#include "idocp/ocp/split_direction.hpp"

namespace idocp {

inline SplitDirection::SplitDirection(const Robot& robot) 
  : split_direction(Eigen::VectorXd::Zero(4*robot.dimv()+robot.dimu())),
    du_passive(Eigen::VectorXd::Zero(robot.dim_passive())),
    dnu_passive(Eigen::VectorXd::Zero(robot.dim_passive())),
    daf_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    dbetamu_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimu_(robot.dimu()), 
    dimx_(2*robot.dimv()), 
    dim_passive_(robot.dim_passive()), 
    dimf_(0), 
    dimKKT_(4*robot.dimv()+robot.dimu()) {
}


inline SplitDirection::SplitDirection() 
  : split_direction(),
    du_passive(),
    dnu_passive(),
    daf_full_(),
    dbetamu_full_(),
    dimv_(0), 
    dimu_(0), 
    dimx_(0), 
    dim_passive_(0), 
    dimf_(0), 
    dimKKT_(0) {
}


inline SplitDirection::~SplitDirection() {
}


inline void SplitDirection::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
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


inline SplitDirection SplitDirection::Random(const Robot& robot) {
  SplitDirection d(robot);
  d.dlmd() = Eigen::VectorXd::Random(robot.dimv());
  d.dgmm() = Eigen::VectorXd::Random(robot.dimv());
  d.du() = Eigen::VectorXd::Random(robot.dimu());
  d.dq() = Eigen::VectorXd::Random(robot.dimv());
  d.dv() = Eigen::VectorXd::Random(robot.dimv());
  d.da() = Eigen::VectorXd::Random(robot.dimv());
  d.dbeta() = Eigen::VectorXd::Random(robot.dimv());
  d.du_passive = Eigen::VectorXd::Random(robot.dim_passive());
  d.dnu_passive = Eigen::VectorXd::Random(robot.dim_passive());
  return d;
}


inline SplitDirection SplitDirection::Random(
    const Robot& robot, const ContactStatus& contact_status) {
  SplitDirection d(robot);
  d.setContactStatus(contact_status);
  d.dlmd() = Eigen::VectorXd::Random(robot.dimv());
  d.dgmm() = Eigen::VectorXd::Random(robot.dimv());
  d.du() = Eigen::VectorXd::Random(robot.dimu());
  d.dq() = Eigen::VectorXd::Random(robot.dimv());
  d.dv() = Eigen::VectorXd::Random(robot.dimv());
  d.da() = Eigen::VectorXd::Random(robot.dimv());
  d.df() = Eigen::VectorXd::Random(contact_status.dimf());
  d.dbeta() = Eigen::VectorXd::Random(robot.dimv());
  d.dmu() = Eigen::VectorXd::Random(contact_status.dimf());
  d.du_passive = Eigen::VectorXd::Random(robot.dim_passive());
  d.dnu_passive = Eigen::VectorXd::Random(robot.dim_passive());
  return d;
}

} // namespace idocp 

#endif // IDOCP_SPLIT_OCP_DIRECTION_HXX_