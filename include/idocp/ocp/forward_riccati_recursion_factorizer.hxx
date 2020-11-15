#ifndef IDOCP_FORWARD_RICCATI_RECURSION_FACTORIZER_HXX_ 
#define IDOCP_FORWARD_RICCATI_RECURSION_FACTORIZER_HXX_

#include "idocp/ocp/forward_riccati_recursion_factorizer.hpp"

#include <cassert>

namespace idocp {

inline ForwardRiccatiRecursionFactorizer::ForwardRiccatiRecursionFactorizer(
    const Robot& robot) 
  : has_floating_base_(robot.has_floating_base()),
    dimv_(robot.dimv()),
    NApBKt_(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())) {
}


inline ForwardRiccatiRecursionFactorizer::ForwardRiccatiRecursionFactorizer() 
  : has_floating_base_(false),
    dimv_(0),
    NApBKt_() {
}


inline ForwardRiccatiRecursionFactorizer::~ForwardRiccatiRecursionFactorizer() {
}


inline void ForwardRiccatiRecursionFactorizer::factorizeStateTransition(
    const RiccatiFactorization& riccati, const KKTMatrix& kkt_matrix, 
    const KKTResidual& kkt_residual, const double dtau,
    RiccatiFactorization& riccati_next) {
  assert(dtau > 0);
  if (has_floating_base_) {
    riccati_next.Pi.topRows(dimv_).noalias()
        = kkt_matrix.Fqq() * riccati.Pi.topRows(dimv_);
    riccati_next.Pi.topRows(dimv_).noalias()
        += dtau * riccati.Pi.bottomRows(dimv_);
  }
  else {
    riccati_next.Pi.topRows(dimv_) 
        = riccati.Pi.topRows(dimv_) + dtau * riccati.Pi.bottomRows(dimv_);
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
  riccati_next.pi.head(dimv_).noalias() += dtau * riccati.pi.tail(dimv_);
  riccati_next.pi.tail(dimv_).noalias() 
      += kkt_matrix.Fxx().bottomRows(dimv_) * riccati.pi;
}


inline void ForwardRiccatiRecursionFactorizer::factorizeStateConstraintFactorization(
    const RiccatiFactorization& riccati, const KKTMatrix& kkt_matrix, 
    const double dtau, RiccatiFactorization& riccati_next) {
  assert(dtau > 0);
  if (has_floating_base_) {
    NApBKt_.leftCols(dimv_).noalias() 
        = riccati.N.leftCols(dimv_) * kkt_matrix.Fqq().transpose();
    NApBKt_.leftCols(dimv_).noalias() 
        += dtau * riccati.N.rightCols(dimv_);
  }
  else {
    NApBKt_.leftCols(dimv_).noalias() 
        = riccati.N.leftCols(dimv_) + dtau * riccati.N.rightCols(dimv_);
  }
  NApBKt_.rightCols(dimv_).noalias() 
      = riccati.N * kkt_matrix.Fxx().bottomRows(dimv_).transpose();
  if (has_floating_base_) {
    riccati_next.N.topRows(dimv_).noalias()
        = kkt_matrix.Fqq() * NApBKt_.topRows(dimv_);
    riccati_next.N.topRows(dimv_).noalias()
        += dtau * NApBKt_.bottomRows(dimv_);
  }
  else {
    riccati_next.N.topRows(dimv_) 
        = NApBKt_.topRows(dimv_) + dtau * NApBKt_.bottomRows(dimv_);
  }
  riccati_next.N.bottomRows(dimv_).noalias() 
      = kkt_matrix.Fxx().bottomRows(dimv_) * NApBKt_;
}

} // namespace idocp

#endif // IDOCP_FORWARD_RICCATI_RECURSION_FACTORIZER_HXX_ 