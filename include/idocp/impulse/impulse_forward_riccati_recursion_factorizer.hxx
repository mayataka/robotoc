#ifndef IDOCP_IMPULSE_FORWARD_RICCATI_RECURSION_FACTORIZER_HXX_ 
#define IDOCP_IMPULSE_FORWARD_RICCATI_RECURSION_FACTORIZER_HXX_

#include "idocp/impulse/impulse_forward_riccati_recursion_factorizer.hpp"

#include <cassert>

namespace idocp {

inline ImpulseForwardRiccatiRecursionFactorizer::ImpulseForwardRiccatiRecursionFactorizer(
    const Robot& robot) 
  : has_floating_base_(robot.has_floating_base()),
    dimv_(robot.dimv()),
    NApBKt_(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())) {
}


inline ImpulseForwardRiccatiRecursionFactorizer::ImpulseForwardRiccatiRecursionFactorizer() 
  : has_floating_base_(false),
    dimv_(0),
    NApBKt_() {
}


inline ImpulseForwardRiccatiRecursionFactorizer::~ImpulseForwardRiccatiRecursionFactorizer() {
}


inline void ImpulseForwardRiccatiRecursionFactorizer::factorizeStateTransition(
    const SplitRiccatiFactorization& riccati, 
    const ImpulseSplitKKTMatrix& kkt_matrix, 
    const ImpulseSplitKKTResidual& kkt_residual,
    SplitRiccatiFactorization& riccati_next) {
  if (has_floating_base_) {
    riccati_next.Pi.topRows(dimv_).noalias()
        = kkt_matrix.Fqq() * riccati.Pi.topRows(dimv_);
  }
  else {
    riccati_next.Pi.topRows(dimv_) = riccati.Pi.topRows(dimv_);
  }
  riccati_next.Pi.bottomRows(dimv_).noalias() 
      = kkt_matrix.Fxx().bottomRows(dimv_) * riccati.Pi;
  riccati_next.pi = kkt_residual.Fx();
  if (has_floating_base_) {
    riccati_next.pi.head(dimv_).noalias() 
        += kkt_matrix.Fqq() * riccati.pi.head(dimv_);
  }
  else {
    riccati_next.pi.head(dimv_).noalias() += riccati.pi.head(dimv_);
  }
  riccati_next.pi.tail(dimv_).noalias() 
      += kkt_matrix.Fxx().bottomRows(dimv_) * riccati.pi;
}


inline void ImpulseForwardRiccatiRecursionFactorizer::factorizeStateConstraintFactorization(
    const SplitRiccatiFactorization& riccati, 
    const ImpulseSplitKKTMatrix& kkt_matrix, 
    SplitRiccatiFactorization& riccati_next) {
  if (has_floating_base_) {
    NApBKt_.leftCols(dimv_).noalias() 
        = riccati.N.leftCols(dimv_) * kkt_matrix.Fqq().transpose();
  }
  else {
    NApBKt_.leftCols(dimv_).noalias() = riccati.N.leftCols(dimv_);
  }
  NApBKt_.rightCols(dimv_).noalias() 
      = riccati.N * kkt_matrix.Fxx().bottomRows(dimv_).transpose();
  if (has_floating_base_) {
    riccati_next.N.topRows(dimv_).noalias()
        = kkt_matrix.Fqq() * NApBKt_.topRows(dimv_);
  }
  else {
    riccati_next.N.topRows(dimv_) = NApBKt_.topRows(dimv_);
  }
  riccati_next.N.bottomRows(dimv_).noalias() 
      = kkt_matrix.Fxx().bottomRows(dimv_) * NApBKt_;
}

} // namespace idocp

#endif // IDOCP_IMPULSE_FORWARD_RICCATI_RECURSION_FACTORIZER_HXX_ 