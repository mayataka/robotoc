#ifndef IDOCP_SPLIT_SWITCHING_CONSTRAINT_RESIDUAL_HXX_ 
#define IDOCP_SPLIT_SWITCHING_CONSTRAINT_RESIDUAL_HXX_

#include "idocp/ocp/split_switching_constraint_residual.hpp"

#include <cassert>


namespace idocp {

inline SplitSwitchingConstraintResidual::SplitSwitchingConstraintResidual(
    const Robot& robot)
  : q(Eigen::VectorXd::Zero(robot.dimq())),
    dq(Eigen::VectorXd::Zero(robot.dimv())),
    P_full_(robot.max_dimf()),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimi_(0) {
}


inline SplitSwitchingConstraintResidual::SplitSwitchingConstraintResidual() 
  : q(),
    dq(),
    P_full_(),
    dimq_(0),
    dimv_(0), 
    dimi_(0) {
}


inline SplitSwitchingConstraintResidual::~SplitSwitchingConstraintResidual() {
}


inline void SplitSwitchingConstraintResidual::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimi_ = impulse_status.dimf();
}


inline void SplitSwitchingConstraintResidual::setImpulseStatus() {
  dimi_ = 0;
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
SplitSwitchingConstraintResidual::P() {
  return P_full_.head(dimi_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitSwitchingConstraintResidual::P() const {
  return P_full_.head(dimi_);
}


inline double SplitSwitchingConstraintResidual::squaredNormKKTResidual() const {
  return P().squaredNorm();
}


inline double SplitSwitchingConstraintResidual::l1NormConstraintViolation() const {
  return P().template lpNorm<1>();
}


inline void SplitSwitchingConstraintResidual::setZero() {
  P().setZero();
}


inline int SplitSwitchingConstraintResidual::dimi() const {
  return dimi_;
}


inline bool SplitSwitchingConstraintResidual::isDimensionConsistent() const {
  if(q.size() != dimq_) return false;
  if(dq.size() != dimv_) return false;
  return true;
}


inline bool SplitSwitchingConstraintResidual::isApprox(
    const SplitSwitchingConstraintResidual& other) const {
  assert(dimi() == other.dimi());
  if (P().isApprox(other.P())) 
    return true;
  else 
    return false;
}


inline bool SplitSwitchingConstraintResidual::hasNaN() const {
  if (P().hasNaN()) 
    return true;
  else 
    return false;
}

} // namespace idocp 

#endif // IDOCP_SPLIT_SWITCHING_CONSTRAINT_RESIDUAL_HXX_ 