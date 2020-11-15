#ifndef IDOCP_IMPULSE_RICCATI_FACTORIZER_HXX_
#define IDOCP_IMPULSE_RICCATI_FACTORIZER_HXX_

#include "idocp/impulse/impulse_riccati_factorizer.hpp"

namespace idocp {

inline ImpulseRiccatiFactorizer::ImpulseRiccatiFactorizer(
    const Robot& robot) 
  : has_floating_base_(robot.has_floating_base()),
    dimv_(robot.dimv()),
    backward_recursion_(robot),
    NApBKt_(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())) {
}


inline ImpulseRiccatiFactorizer::ImpulseRiccatiFactorizer() 
  : has_floating_base_(false),
    dimv_(0),
    backward_recursion_(),
    NApBKt_() {
}


inline ImpulseRiccatiFactorizer::~ImpulseRiccatiFactorizer() {
}


inline void ImpulseRiccatiFactorizer::backwardRiccatiRecursion(
    const RiccatiFactorization& riccati_next, ImpulseKKTMatrix& kkt_matrix, 
    ImpulseKKTResidual& kkt_residual, RiccatiFactorization& riccati) {
  backward_recursion_.factorizeKKTMatrix(riccati_next, kkt_matrix, kkt_residual);
  backward_recursion_.factorizeRiccatiFactorization(riccati_next, kkt_matrix, 
                                                    kkt_residual, riccati);
}


inline void ImpulseRiccatiFactorizer::forwardRiccatiRecursionSerial(
    const RiccatiFactorization& riccati, const ImpulseKKTMatrix& kkt_matrix, 
    const ImpulseKKTResidual& kkt_residual, 
    RiccatiFactorization& riccati_next) {
  riccati_next.Pi.noalias() = kkt_matrix.Fxx() * riccati.Pi;
  riccati_next.pi = kkt_residual.Fx();
  riccati_next.pi.noalias() += kkt_matrix.Fxx() * riccati.pi;
  NApBKt_.noalias() = riccati.N * kkt_matrix.Fxx().transpose();
  riccati_next.N.noalias() = kkt_matrix.Fxx() * NApBKt_;
}


template <typename VectorType>
inline void ImpulseRiccatiFactorizer::computeStateDirection(
    const RiccatiFactorization& riccati, 
    const Eigen::MatrixBase<VectorType>& dx0, ImpulseSplitDirection& d) {
  d.dx().noalias() = riccati.Pi * dx0;
  d.dx().noalias() += riccati.pi;
  d.dx().noalias() -= riccati.N * riccati.n;
}


inline void ImpulseRiccatiFactorizer::computeCostateDirection(
    const RiccatiFactorization& riccati, ImpulseSplitDirection& d) {
  d.dlmd().noalias() = riccati.Pqq * d.dq();
  d.dlmd().noalias() += riccati.Pqv * d.dv();
  d.dlmd().noalias() -= riccati.sq;
  d.dgmm().noalias() = riccati.Pqv.transpose() * d.dq();
  d.dgmm().noalias() += riccati.Pvv * d.dv();
  d.dgmm().noalias() -= riccati.sv;
  d.dlmdgmm().noalias() += riccati.n;
}

} // namespace idocp

#endif // IDOCP_IMPULSE_RICCATI_FACTORIZER_HXX_