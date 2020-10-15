#ifndef IDOCP_IMPULSE_SPLIT_DIRECTION_HXX_
#define IDOCP_IMPULSE_SPLIT_DIRECTION_HXX_

#include "idocp/impulse/impulse_split_direction.hpp"

namespace idocp {

inline ImpulseSplitDirection::ImpulseSplitDirection(
    const Robot& robot, const bool use_contact_position_constraint) 
  : ddv(robot.dimv()),
    dbeta(robot.dimv()),
    split_direction_(Eigen::VectorXd::Zero(4*robot.dimv()+3*robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimf_(0), 
    dimc_(0),
    max_dimKKT_(4*robot.dimv()+3*robot.max_dimf()),
    dimKKT_(4*robot.dimv()) {
}


inline ImpulseSplitDirection::ImpulseSplitDirection() 
  : ddv(),
    dbeta(),
    split_direction_(),
    dimv_(0), 
    dimx_(0), 
    dimf_(0), 
    dimc_(0),
    max_dimKKT_(0),
    dimKKT_(0) {
}


inline ImpulseSplitDirection::~ImpulseSplitDirection() {
}


inline void ImpulseSplitDirection::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
  dimc_ = 2*contact_status.dimf();
  dimKKT_ = 4*dimv_ + dimf_ + dimc_;
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseSplitDirection::split_direction() {
  return split_direction_.head(dimKKT_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dlmd() {
  return split_direction_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dgmm() {
  return split_direction_.segment(dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dmu() {
  return split_direction_.segment(dimx_, dimc_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseSplitDirection::dmu_contact_position() {
  return split_direction_.segment(dimx_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseSplitDirection::dmu_contact_velocity() {
  return split_direction_.segment(dimx_+dimf_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::df() {
  return split_direction_.segment(dimx_+dimc_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dq() {
  return split_direction_.segment(dimx_+dimc_+dimf_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dv() {
  return split_direction_.segment(dimx_+dimc_+dimf_+dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitDirection::dx() {
  return split_direction_.segment(dimx_+dimc_+dimf_, dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::split_direction() const {
  return split_direction_.head(dimKKT_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dlmd() const {
  return split_direction_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dgmm() const {
  return split_direction_.segment(dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dmu() const {
  return split_direction_.segment(dimx_, dimc_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dmu_contact_position() const {
  return split_direction_.segment(dimx_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dmu_contact_velocity() const {
  return split_direction_.segment(dimx_+dimf_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::df() const {
  return split_direction_.segment(dimx_+dimc_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dq() const {
  return split_direction_.segment(dimx_+dimc_+dimf_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dv() const {
  return split_direction_.segment(dimx_+dimc_+dimf_+dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitDirection::dx() const {
  return split_direction_.segment(dimx_+dimc_+dimf_, dimx_);
}


inline void ImpulseSplitDirection::setZero() {
  ddv.setZero();
  dbeta.setZero();
  split_direction_.setZero();
}


inline int ImpulseSplitDirection::dimKKT() const {
  return dimKKT_;
}


inline int ImpulseSplitDirection::max_dimKKT() const {
  return max_dimKKT_;
}


inline int ImpulseSplitDirection::dimc() const {
  return dimc_;
}


inline int ImpulseSplitDirection::dimf() const {
  return dimf_;
}


inline ImpulseSplitDirection ImpulseSplitDirection::Random(
    const Robot& robot, const ContactStatus& contact_status, 
    const bool use_contact_position_constraint) {
  ImpulseSplitDirection d(robot, use_contact_position_constraint);
  d.setContactStatus(contact_status);
  d.dlmd() = Eigen::VectorXd::Random(robot.dimv());
  d.dgmm() = Eigen::VectorXd::Random(robot.dimv());
  d.dmu() = Eigen::VectorXd::Random(d.dimc());
  d.df() = Eigen::VectorXd::Random(d.dimf());
  d.dq() = Eigen::VectorXd::Random(robot.dimv());
  d.dv() = Eigen::VectorXd::Random(robot.dimv());
  d.ddv = Eigen::VectorXd::Random(robot.dimv());
  d.dbeta = Eigen::VectorXd::Random(robot.dimv());
  return d;
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_SPLIT_DIRECTION_HXX_ 