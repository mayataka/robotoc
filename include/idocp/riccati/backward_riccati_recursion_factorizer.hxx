#ifndef IDOCP_BACKWARD_RICCATI_RECURSION_FACTORIZER_HXX_ 
#define IDOCP_BACKWARD_RICCATI_RECURSION_FACTORIZER_HXX_

#include "idocp/riccati/backward_riccati_recursion_factorizer.hpp"

namespace idocp {

inline BackwardRiccatiRecursionFactorizer::BackwardRiccatiRecursionFactorizer(
    const Robot& robot) 
  : dimv_(robot.dimv()),
    dimu_(robot.dimu()),
    AtP_(MatrixXdRowMajor::Zero(2*robot.dimv(), 2*robot.dimv())),
    BtP_(MatrixXdRowMajor::Zero(robot.dimu(), 2*robot.dimv())),
    GK_(Eigen::MatrixXd::Zero(robot.dimu(), 2*robot.dimv())) {
}


inline BackwardRiccatiRecursionFactorizer::BackwardRiccatiRecursionFactorizer() 
  : dimv_(0),
    dimu_(0),
    AtP_(),
    BtP_(),
    GK_() {
}


inline BackwardRiccatiRecursionFactorizer::
~BackwardRiccatiRecursionFactorizer() {
}


inline void BackwardRiccatiRecursionFactorizer::factorizeKKTMatrix(
    const SplitRiccatiFactorization& riccati_next, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) {
  AtP_.noalias() = kkt_matrix.Fxx.transpose() * riccati_next.P;
  BtP_.noalias() = kkt_matrix.Fvu.transpose() * riccati_next.P.bottomRows(dimv_);
  // Factorize F
  kkt_matrix.Qxx.noalias() += AtP_ * kkt_matrix.Fxx;
  // Factorize H
  kkt_matrix.Qxu.noalias() += AtP_.rightCols(dimv_) * kkt_matrix.Fvu;
  // Factorize G
  kkt_matrix.Quu.noalias() += BtP_.rightCols(dimv_) * kkt_matrix.Fvu;
  // Factorize vector term
  kkt_residual.lu.noalias() += BtP_ * kkt_residual.Fx;
  kkt_residual.lu.noalias() -= kkt_matrix.Fvu.transpose() * riccati_next.sv();
}


inline void BackwardRiccatiRecursionFactorizer::factorizeKKTMatrix(
    const SplitRiccatiFactorization& riccati_next, 
    ImpulseSplitKKTMatrix& kkt_matrix) {
  AtP_.noalias() = kkt_matrix.Fxx.transpose() * riccati_next.P;
  // Factorize F
  kkt_matrix.Qxx.noalias() += AtP_ * kkt_matrix.Fxx;
}


inline void BackwardRiccatiRecursionFactorizer::factorizeRiccatiFactorization(
    const SplitRiccatiFactorization& riccati_next, SplitKKTMatrix& kkt_matrix, 
    const SplitKKTResidual& kkt_residual, const LQRPolicy& lqr_policy, 
    SplitRiccatiFactorization& riccati) {
  GK_.noalias() = kkt_matrix.Quu * lqr_policy.K; 
  kkt_matrix.Qxx.noalias() -= lqr_policy.K.transpose() * GK_;
  // Riccati factorization matrix with preserving the symmetry
  riccati.P = 0.5 * (kkt_matrix.Qxx + kkt_matrix.Qxx.transpose());
  // Riccati factorization vector
  riccati.s.noalias()  = kkt_matrix.Fxx.transpose() * riccati_next.s;
  riccati.s.noalias() -= AtP_ * kkt_residual.Fx;
  riccati.s.noalias() -= kkt_residual.lx;
  riccati.s.noalias() -= kkt_matrix.Qxu * lqr_policy.k;
}


inline void BackwardRiccatiRecursionFactorizer::factorizeRiccatiFactorization(
    const SplitRiccatiFactorization& riccati_next, 
    const ImpulseSplitKKTMatrix& kkt_matrix, 
    const ImpulseSplitKKTResidual& kkt_residual, 
    SplitRiccatiFactorization& riccati) {
  // Riccati factorization matrix with preserving the symmetry
  riccati.P = 0.5 * (kkt_matrix.Qxx + kkt_matrix.Qxx.transpose());
  // Riccati factorization vector
  riccati.s.noalias()  = kkt_matrix.Fxx.transpose() * riccati_next.s;
  riccati.s.noalias() -= AtP_ * kkt_residual.Fx;
  riccati.s.noalias() -= kkt_residual.lx;
}

} // namespace idocp

#endif // IDOCP_BACKWARD_RICCATI_RECURSION_FACTORIZER_HXX_ 