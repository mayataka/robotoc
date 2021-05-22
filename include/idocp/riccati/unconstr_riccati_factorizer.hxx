#ifndef IDOCP_UNCONSTR_RICCATI_FACTORIZER_HXX_ 
#define IDOCP_UNCONSTR_RICCATI_FACTORIZER_HXX_

#include "idocp/riccati/unconstr_riccati_factorizer.hpp"

#include <cassert>

namespace idocp {

inline UnconstrRiccatiFactorizer::UnconstrRiccatiFactorizer(const Robot& robot) 
  : dimv_(robot.dimv()),
    llt_(robot.dimv()),
    backward_recursion_(robot) {
}


inline UnconstrRiccatiFactorizer::UnconstrRiccatiFactorizer() 
  : dimv_(0),
    llt_(),
    backward_recursion_() {
}


inline UnconstrRiccatiFactorizer::~UnconstrRiccatiFactorizer() {
}


inline void UnconstrRiccatiFactorizer::backwardRiccatiRecursion(
    const SplitRiccatiFactorization& riccati_next, const double dt, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual,  
    SplitRiccatiFactorization& riccati, LQRPolicy& lqr_policy) {
  assert(dt > 0);
  backward_recursion_.factorizeKKTMatrix(riccati_next, dt, kkt_matrix, kkt_residual);
  llt_.compute(kkt_matrix.Qaa);
  assert(llt_.info() == Eigen::Success);
  lqr_policy.K = - llt_.solve(kkt_matrix.Qxu.transpose());
  lqr_policy.k = - llt_.solve(kkt_residual.la);
  assert(!lqr_policy.K.hasNaN());
  assert(!lqr_policy.k.hasNaN());
  backward_recursion_.factorizeRiccatiFactorization(riccati_next, kkt_matrix, 
                                                    kkt_residual, lqr_policy,
                                                    dt, riccati);
}


inline void UnconstrRiccatiFactorizer::forwardRiccatiRecursion(
    const SplitKKTResidual& kkt_residual, const double dt, 
    const LQRPolicy& lqr_policy, SplitDirection& d, 
    SplitDirection& d_next) const {
  assert(dt > 0);
  d.da().noalias() = lqr_policy.K * d.dx + lqr_policy.k;
  d_next.dx = kkt_residual.Fx + d.dx;
  d_next.dq().noalias() += dt * d.dv();
  d_next.dv().noalias() += dt * d.da();
}


inline void UnconstrRiccatiFactorizer::computeCostateDirection(
    const SplitRiccatiFactorization& riccati, SplitDirection& d) {
  d.dlmdgmm.noalias() = riccati.P * d.dx - riccati.s;
}

} // namespace idocp

#endif // IDOCP_UNCONSTR_RICCATI_FACTORIZER_HXX_ 