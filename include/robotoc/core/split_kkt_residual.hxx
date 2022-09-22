#ifndef ROBOTOC_SPLIT_KKT_RESIDUAL_HXX_ 
#define ROBOTOC_SPLIT_KKT_RESIDUAL_HXX_

#include "robotoc/core/split_kkt_residual.hpp"

#include <cmath>


namespace robotoc {

inline void SplitKKTResidual::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
}


inline void SplitKKTResidual::setContactStatus(
    const ImpulseStatus& contact_status) {
  dimf_ = contact_status.dimf();
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
  err += lx.squaredNorm();
  err += lu.squaredNorm();
  err += la.squaredNorm();
  err += ldv.squaredNorm();
  err += lf().squaredNorm();
  return err;
}


template <int p>
inline double SplitKKTResidual::constraintViolation() const {
  return Fx.template lpNorm<p>();
}


inline void SplitKKTResidual::setZero() {
  Fx.setZero();
  lx.setZero();
  la.setZero();
  ldv.setZero();
  lu.setZero();
  lf().setZero();
  h = 0.0;
  kkt_error = 0.0;
  cost = 0.0;
  constraint_violation = 0.0;
}


inline int SplitKKTResidual::dimf() const {
  return dimf_;
}

} // namespace robotoc 

#endif // ROBOTOC_SPLIT_KKT_RESIDUAL_HXX_ 