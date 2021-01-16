#ifndef IDOCP_PRE_IMPULSE_SPLIT_KKT_RESIDUAL_HXX_ 
#define IDOCP_PRE_IMPULSE_SPLIT_KKT_RESIDUAL_HXX_

#include "idocp/ocp/pre_impulse_split_kkt_residual.hpp"


namespace idocp {

inline PreImpulseSplitKKTResidual::PreImpulseSplitKKTResidual(const Robot& robot) 
  : P_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    dimf_(0) {
}


inline PreImpulseSplitKKTResidual::PreImpulseSplitKKTResidual() 
  : P_full_(),
    dimf_(0) {
}


inline PreImpulseSplitKKTResidual::~PreImpulseSplitKKTResidual() {
}


inline void PreImpulseSplitKKTResidual::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimf_ = impulse_status.dimf();
}


inline Eigen::VectorBlock<Eigen::VectorXd> PreImpulseSplitKKTResidual::P() {
  return P_full_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
PreImpulseSplitKKTResidual::P() const {
  return P_full_.head(dimf_);
}


inline void PreImpulseSplitKKTResidual::setZero() {
  P_full_.setZero();
}


inline int PreImpulseSplitKKTResidual::dimf() const {
  return dimf_;
}


inline bool PreImpulseSplitKKTResidual::isApprox(
    const PreImpulseSplitKKTResidual& other) const {
  if (!P().isApprox(other.P())) return false;
  return true;
}


inline bool PreImpulseSplitKKTResidual::hasNaN() const {
  if (P_full_.hasNaN()) return true;
  return false;
}

} // namespace idocp 

#endif // IDOCP_PRE_IMPULSE_SPLIT_KKT_RESIDUAL_HXX_ 