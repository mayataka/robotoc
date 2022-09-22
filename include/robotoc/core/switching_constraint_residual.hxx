#ifndef ROBOTOC_SWITCHING_CONSTRAINT_RESIDUAL_HXX_ 
#define ROBOTOC_SWITCHING_CONSTRAINT_RESIDUAL_HXX_

#include "robotoc/core/switching_constraint_residual.hpp"

#include <cassert>


namespace robotoc {

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

} // namespace robotoc 

#endif // ROBOTOC_SWITCHING_CONSTRAINT_RESIDUAL_HXX_ 