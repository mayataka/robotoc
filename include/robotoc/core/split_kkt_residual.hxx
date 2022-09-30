#ifndef ROBOTOC_SPLIT_KKT_RESIDUAL_HXX_ 
#define ROBOTOC_SPLIT_KKT_RESIDUAL_HXX_

#include "robotoc/core/split_kkt_residual.hpp"

#include <cmath>


namespace robotoc {

inline void SplitKKTResidual::setContactDimension(const int dimf) {
  assert(dimf >= 0);
  assert(dimf <= P_full_.size());
  dimf_ = dimf;
}


inline void SplitKKTResidual::setSwitchingConstraintDimension(const int dims) {
  assert(dims >= 0);
  assert(dims <= P_full_.size());
  dims_ = dims;
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::Fq() {
  return Fx.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::Fq() const {
  return Fx.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::Fv() {
  return Fx.tail(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::Fv() const {
  return Fx.tail(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::P() {
  return P_full_.head(dims_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> SplitKKTResidual::P() const {
  return P_full_.head(dims_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::lq() {
  return lx.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::lq() const {
  return lx.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::lv() {
  return lx.tail(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::lv() const {
  return lx.tail(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::lf() {
  return lf_full_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::lf() const {
  return lf_full_.head(dimf_);
}


inline double SplitKKTResidual::KKTError() const {
  double err = 0;
  err += Fx.squaredNorm();
  if (P().size() > 0) {
    err += P().squaredNorm();
  }
  err += lx.squaredNorm();
  err += lu.squaredNorm();
  err += la.squaredNorm();
  err += ldv.squaredNorm();
  if (lf().size() > 0) {
    err += lf().squaredNorm();
  }
  return err;
}


template <int p>
inline double SplitKKTResidual::constraintViolation() const {
  double vio = Fx.template lpNorm<p>();
  if (P().size() > 0) {
    vio += P().template lpNorm<p>();
  }
  return vio;
}


inline void SplitKKTResidual::setZero() {
  Fx.setZero();
  if (P().size() > 0) {
    P().setZero();
  }
  lx.setZero();
  la.setZero();
  ldv.setZero();
  lu.setZero();
  if (lf().size() > 0) {
    lf().setZero();
  }
  h = 0.0;
  kkt_error = 0.0;
  cost = 0.0;
  constraint_violation = 0.0;
}


inline int SplitKKTResidual::dimf() const {
  return dimf_;
}


inline int SplitKKTResidual::dims() const {
  return dims_;
}

} // namespace robotoc 

#endif // ROBOTOC_SPLIT_KKT_RESIDUAL_HXX_ 