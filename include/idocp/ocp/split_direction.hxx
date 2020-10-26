#ifndef IDOCP_SPLIT_DIRECTION_HXX_
#define IDOCP_SPLIT_DIRECTION_HXX_

#include "idocp/ocp/split_direction.hpp"

namespace idocp {

inline SplitDirection::SplitDirection(const Robot& robot) 
  : split_direction(Eigen::VectorXd::Zero(4*robot.dimv()+robot.dimu())),
    da(Eigen::VectorXd::Zero(robot.dimv())),
    dbeta(Eigen::VectorXd::Zero(robot.dimv())),
    dnu_passive(Eigen::VectorXd::Zero(robot.dim_passive())),
    dmu_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    df_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimu_(robot.dimu()), 
    dimx_(2*robot.dimv()), 
    dim_passive_(robot.dim_passive()), 
    dimf_(0), 
    dimKKT_(4*robot.dimv()+robot.dimu()) {
}


inline SplitDirection::SplitDirection() 
  : split_direction(),
    da(),
    dbeta(),
    dnu_passive(),
    dmu_full_(),
    df_full_(),
    dimv_(), 
    dimu_(), 
    dimx_(), 
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


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dmu() {
  return dmu_full_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dmu() const {
  return dmu_full_.head(dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::df() {
  return df_full_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::df() const {
  return df_full_.head(dimf_);
}


inline void SplitDirection::setZero() {
  split_direction.setZero();
  da.setZero();
  dbeta.setZero();
  du_passive.setZero();
  dnu_passive.setZero();
  dmu_full_.setZero();
  df_full_.setZero();
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
  d.da = Eigen::VectorXd::Random(robot.dimv());
  d.dbeta = Eigen::VectorXd::Random(robot.dimv());
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
  d.dmu() = Eigen::VectorXd::Random(contact_status.dimf());
  d.df() = Eigen::VectorXd::Random(contact_status.dimf());
  d.da = Eigen::VectorXd::Random(robot.dimv());
  d.dbeta = Eigen::VectorXd::Random(robot.dimv());
  d.du_passive = Eigen::VectorXd::Random(robot.dim_passive());
  d.dnu_passive = Eigen::VectorXd::Random(robot.dim_passive());
  return d;
}

} // namespace idocp 

#endif // IDOCP_SPLIT_OCP_DIRECTION_HXX_