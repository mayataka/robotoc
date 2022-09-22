#ifndef ROBOTOC_SWITCHING_CONSTRAINT_RESIDUAL_HXX_ 
#define ROBOTOC_SWITCHING_CONSTRAINT_RESIDUAL_HXX_

#include "robotoc/core/switching_constraint_residual.hpp"

#include <cassert>


namespace robotoc {

inline void SwitchingConstraintResidual::setDimension(const int dims) {
  assert(dims >= 0);
  assert(dims <= P_full_.size());
  dims_ = dims;
}


inline Eigen::VectorBlock<Eigen::VectorXd> SwitchingConstraintResidual::P() {
  return P_full_.head(dims_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SwitchingConstraintResidual::P() const {
  return P_full_.head(dims_);
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


inline int SwitchingConstraintResidual::dims() const {
  return dims_;
}

} // namespace robotoc 

#endif // ROBOTOC_SWITCHING_CONSTRAINT_RESIDUAL_HXX_ 