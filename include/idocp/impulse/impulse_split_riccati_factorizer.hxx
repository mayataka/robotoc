#ifndef IDOCP_IMPULSE_SPLIT_RICCATI_FACTORIZER_HXX_
#define IDOCP_IMPULSE_SPLIT_RICCATI_FACTORIZER_HXX_

#include "idocp/impulse/impulse_split_riccati_factorizer.hpp"

namespace idocp {

inline ImpulseSplitRiccatiFactorizer::ImpulseSplitRiccatiFactorizer(
    const Robot& robot) 
  : has_floating_base_(robot.hasFloatingBase()),
    dimv_(robot.dimv()),
    backward_recursion_(robot) {
}


inline ImpulseSplitRiccatiFactorizer::ImpulseSplitRiccatiFactorizer() 
  : has_floating_base_(false),
    dimv_(0),
    backward_recursion_() {
}


inline ImpulseSplitRiccatiFactorizer::~ImpulseSplitRiccatiFactorizer() {
}


inline void ImpulseSplitRiccatiFactorizer::backwardRiccatiRecursion(
    const SplitRiccatiFactorization& riccati_next, 
    ImpulseSplitKKTMatrix& kkt_matrix, ImpulseSplitKKTResidual& kkt_residual, 
    SplitRiccatiFactorization& riccati) {
  backward_recursion_.factorizeKKTMatrix(riccati_next, kkt_matrix);
  backward_recursion_.factorizeRiccatiFactorization(riccati_next, kkt_matrix, 
                                                    kkt_residual, riccati);
}


inline void ImpulseSplitRiccatiFactorizer::forwardRiccatiRecursion(
    const ImpulseSplitKKTMatrix& kkt_matrix, 
    const ImpulseSplitKKTResidual& kkt_residual, 
    const ImpulseSplitDirection& d, SplitDirection& d_next) const {
  d_next.dx = kkt_residual.Fx;
  if (has_floating_base_) {
    d_next.dq().template head<6>().noalias() 
        += kkt_matrix.Fqq().template topLeftCorner<6, 6>() 
            * d.dq().template head<6>();
    d_next.dq().tail(dimv_-6).noalias() += d.dq().tail(dimv_-6);
  }
  else {
    d_next.dq().noalias() += d.dq();
  }
  d_next.dv().noalias() += kkt_matrix.Fvq() * d.dq();
  d_next.dv().noalias() += kkt_matrix.Fvv() * d.dv();
}


inline void ImpulseSplitRiccatiFactorizer::computeCostateDirection(
    const SplitRiccatiFactorization& riccati, ImpulseSplitDirection& d) {
  d.dlmd().noalias()  = riccati.Pqq * d.dq();
  d.dlmd().noalias() += riccati.Pqv * d.dv();
  d.dlmd().noalias() -= riccati.sq;
  d.dgmm().noalias()  = riccati.Pqv.transpose() * d.dq();
  d.dgmm().noalias() += riccati.Pvv * d.dv();
  d.dgmm().noalias() -= riccati.sv;
}

} // namespace idocp

#endif // IDOCP_IMPULSE_SPLIT_RICCATI_FACTORIZER_HXX_ 