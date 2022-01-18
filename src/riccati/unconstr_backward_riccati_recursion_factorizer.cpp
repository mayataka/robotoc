#include "robotoc/riccati/unconstr_backward_riccati_recursion_factorizer.hpp"

#include <cassert>


namespace robotoc {

UnconstrBackwardRiccatiRecursionFactorizer::
UnconstrBackwardRiccatiRecursionFactorizer(const Robot& robot) 
  : dimv_(robot.dimv()),
    GK_(Eigen::MatrixXd::Zero(robot.dimv(), 2*robot.dimv())) {
}


UnconstrBackwardRiccatiRecursionFactorizer::
UnconstrBackwardRiccatiRecursionFactorizer() 
  : dimv_(0),
    GK_() {
}


UnconstrBackwardRiccatiRecursionFactorizer::
~UnconstrBackwardRiccatiRecursionFactorizer() {
}


void UnconstrBackwardRiccatiRecursionFactorizer::factorizeKKTMatrix(
    const SplitRiccatiFactorization& riccati_next, const double dt, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  // Factorize F
  kkt_matrix.Qxx.noalias() += riccati_next.P;
  kkt_matrix.Qxx.bottomRows(dimv_).noalias() 
      += dt * riccati_next.P.topRows(dimv_);
  kkt_matrix.Qxx.rightCols(dimv_).noalias() 
      += dt * riccati_next.P.leftCols(dimv_);
  kkt_matrix.Qxx.bottomRightCorner(dimv_, dimv_).noalias() 
      += (dt*dt) * riccati_next.P.topLeftCorner(dimv_, dimv_);
  // Factorize H
  kkt_matrix.Qxu.noalias()   += dt * riccati_next.P.rightCols(dimv_);
  kkt_matrix.Qvu().noalias() += (dt*dt) * riccati_next.Pqv();
  // Factorize G
  kkt_matrix.Qaa.noalias()   += (dt*dt) * riccati_next.Pvv();
  // Factorize vector term
  kkt_residual.la.noalias() += dt * riccati_next.P.bottomRows(dimv_) 
                                  * kkt_residual.Fx;
  kkt_residual.la.noalias() -= dt * riccati_next.sv();
}


void UnconstrBackwardRiccatiRecursionFactorizer::
factorizeRiccatiFactorization(const SplitRiccatiFactorization& riccati_next, 
                              SplitKKTMatrix& kkt_matrix, 
                              const SplitKKTResidual& kkt_residual, 
                              const LQRPolicy& lqr_policy, const double dt, 
                              SplitRiccatiFactorization& riccati) {
  assert(dt > 0);
  GK_.noalias() = kkt_matrix.Qaa * lqr_policy.K; 
  kkt_matrix.Qxx.noalias() -= lqr_policy.K.transpose() * GK_;
  // Riccati factorization matrix with preserving the symmetry
  riccati.P = 0.5 * (kkt_matrix.Qxx + kkt_matrix.Qxx.transpose());
  // Riccati factorization vector
  riccati.s = riccati_next.s;
  riccati.sv().noalias() += dt * riccati_next.sq();
  riccati.s.noalias() -= riccati_next.P * kkt_residual.Fx;
  riccati.sv().noalias() -= dt * riccati_next.P.topRows(dimv_) 
                               * kkt_residual.Fx;
  riccati.s.noalias() -= kkt_residual.lx;
  riccati.s.noalias() -= kkt_matrix.Qxu * lqr_policy.k;
}

} // namespace robotoc
