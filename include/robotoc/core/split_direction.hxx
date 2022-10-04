#ifndef ROBOTOC_SPLIT_DIRECTION_HXX_
#define ROBOTOC_SPLIT_DIRECTION_HXX_

#include "robotoc/core/split_direction.hpp"

#include <cassert>

namespace robotoc {

inline void SplitDirection::setContactDimension(const int dimf) {
  assert(dimf >= 0);
  assert(dimf <= dxi_full_.size());
  dimf_ = dimf;
}


inline void SplitDirection::setSwitchingConstraintDimension(const int dims) {
  assert(dims >= 0);
  assert(dims <= dxi_full_.size());
  dims_ = dims;
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
  return dxi_full_.head(dims_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dxi() const {
  return dxi_full_.head(dims_);
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


inline int SplitDirection::dims() const {
  return dims_;
}

} // namespace robotoc 

#endif // ROBOTOC_SPLIT_OCP_DIRECTION_HXX_