#ifndef ROBOTOC_SWITCHING_CONSTRAINT_RESIDUAL_HXX_ 
#define ROBOTOC_SWITCHING_CONSTRAINT_RESIDUAL_HXX_

#include "robotoc/ocp/switching_constraint_residual.hpp"

#include <cassert>


namespace robotoc {

inline SwitchingConstraintResidual::SwitchingConstraintResidual(
    const Robot& robot)
  : P_full_(robot.max_dimf()),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimi_(0) {
}


inline SwitchingConstraintResidual::SwitchingConstraintResidual() 
  : P_full_(),
    dimq_(0),
    dimv_(0), 
    dimi_(0) {
}


inline SwitchingConstraintResidual::~SwitchingConstraintResidual() {
}


inline void SwitchingConstraintResidual::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimi_ = impulse_status.dimi();
}


inline void SwitchingConstraintResidual::setImpulseStatus() {
  dimi_ = 0;
}


inline Eigen::VectorBlock<Eigen::VectorXd> SwitchingConstraintResidual::P() {
  return P_full_.head(dimi_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SwitchingConstraintResidual::P() const {
  return P_full_.head(dimi_);
}


inline double SwitchingConstraintResidual::KKTError() const {
  return P().squaredNorm();
}


template <int p>
inline double SwitchingConstraintResidual::constraintViolation() const {
  return P().template lpNorm<p>();
}


inline void SwitchingConstraintResidual::setZero() {
  P().setZero();
}


inline int SwitchingConstraintResidual::dimi() const {
  return dimi_;
}


inline bool SwitchingConstraintResidual::isApprox(
    const SwitchingConstraintResidual& other) const {
  assert(dimi() == other.dimi());
  if (P().isApprox(other.P())) 
    return true;
  else 
    return false;
}


inline bool SwitchingConstraintResidual::hasNaN() const {
  if (P().hasNaN()) 
    return true;
  else 
    return false;
}

} // namespace robotoc 

#endif // ROBOTOC_SWITCHING_CONSTRAINT_RESIDUAL_HXX_ 