#ifndef IDOCP_IMPULSE_FORWARD_RICCATI_RECURSION_FACTORIZER_HXX_ 
#define IDOCP_IMPULSE_FORWARD_RICCATI_RECURSION_FACTORIZER_HXX_

#include "idocp/impulse/impulse_forward_riccati_recursion_factorizer.hpp"

#include <cassert>

namespace idocp {

inline ImpulseForwardRiccatiRecursionFactorizer::
ImpulseForwardRiccatiRecursionFactorizer(const Robot& robot) 
  : has_floating_base_(robot.hasFloatingBase()),
    dimv_(robot.dimv()),
    NApBKt_(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())) {
}


inline ImpulseForwardRiccatiRecursionFactorizer::
ImpulseForwardRiccatiRecursionFactorizer() 
  : has_floating_base_(false),
    dimv_(0),
    NApBKt_() {
}


inline ImpulseForwardRiccatiRecursionFactorizer::
~ImpulseForwardRiccatiRecursionFactorizer() {
}


inline void ImpulseForwardRiccatiRecursionFactorizer::factorizeStateTransition(
    const SplitRiccatiFactorization& riccati, 
    const ImpulseSplitKKTMatrix& kkt_matrix, 
    const ImpulseSplitKKTResidual& kkt_residual,
    SplitRiccatiFactorization& riccati_next) {
  // if (has_floating_base_) {
  //   riccati_next.Pi.template topRows<6>().noalias() 
  //       = kkt_matrix.Fqq().template topLeftCorner<6, 6>() 
  //           * riccati.Pi.template topRows<6>();
  //   riccati_next.Pi.middleRows(6, dimv_-6) = riccati.Pi.middleRows(6, dimv_-6);
  // }
  // else {
  //   riccati_next.Pi.topRows(dimv_) = riccati.Pi.topRows(dimv_);
  // }
  // riccati_next.Pi.bottomRows(dimv_).noalias() 
  //     = kkt_matrix.Fxx().bottomRows(dimv_) * riccati.Pi;
  // riccati_next.pi = kkt_residual.Fx();
  // if (has_floating_base_) {
  //   riccati_next.pi.template head<6>().noalias() 
  //       += kkt_matrix.Fqq().template topLeftCorner<6, 6>() 
  //           * riccati.pi.template head<6>();
  //   riccati_next.pi.segment(6, dimv_-6).noalias() 
  //       += riccati.pi.segment(6, dimv_-6);
  // }
  // else {
  //   riccati_next.pi.head(dimv_).noalias() += riccati.pi.head(dimv_);
  // }
  // riccati_next.pi.tail(dimv_).noalias() 
  //     += kkt_matrix.Fxx().bottomRows(dimv_) * riccati.pi;
}


inline void ImpulseForwardRiccatiRecursionFactorizer::
factorizeStateConstraintFactorization(
    const SplitRiccatiFactorization& riccati, 
    const ImpulseSplitKKTMatrix& kkt_matrix, 
    SplitRiccatiFactorization& riccati_next) {
  // if (has_floating_base_) {
  //   NApBKt_.template leftCols<6>().noalias() 
  //       = riccati.N.template leftCols<6>() 
  //           * kkt_matrix.Fqq().template topLeftCorner<6, 6>().transpose();
  //   NApBKt_.middleCols(6, dimv_-6) = riccati.N.middleCols(6, dimv_-6);
  // }
  // else {
  //   NApBKt_.leftCols(dimv_) = riccati.N.leftCols(dimv_);
  // }
  // NApBKt_.rightCols(dimv_).noalias() 
  //     = riccati.N * kkt_matrix.Fxx().bottomRows(dimv_).transpose();
  // if (has_floating_base_) {
  //   riccati_next.N.template topRows<6>().noalias() 
  //       = kkt_matrix.Fqq().template topLeftCorner<6, 6>() 
  //           * NApBKt_.template topRows<6>();
  //   riccati_next.N.middleRows(6, dimv_-6) = NApBKt_.middleRows(6, dimv_-6);
  // }
  // else {
  //   riccati_next.N.topRows(dimv_) = NApBKt_.topRows(dimv_);
  // }
  // riccati_next.N.bottomRows(dimv_).noalias() 
  //     = kkt_matrix.Fxx().bottomRows(dimv_) * NApBKt_;
}

} // namespace idocp

#endif // IDOCP_IMPULSE_FORWARD_RICCATI_RECURSION_FACTORIZER_HXX_ 