#ifndef IDOCP_BACKWARD_UNRICCATI_RECURSION_FACTORIZER_HXX_ 
#define IDOCP_BACKWARD_UNRICCATI_RECURSION_FACTORIZER_HXX_

#include "idocp/unocp/backward_unriccati_recursion_factorizer.hpp"

#include <cassert>

namespace idocp {

inline BackwardUnRiccatiRecursionFactorizer::
BackwardUnRiccatiRecursionFactorizer(const Robot& robot) 
  : dimv_(robot.dimv()),
    GK_(Eigen::MatrixXd::Zero(robot.dimv(), 2*robot.dimv())) {
}


inline BackwardUnRiccatiRecursionFactorizer::
BackwardUnRiccatiRecursionFactorizer() 
  : dimv_(0),
    GK_() {
}


inline BackwardUnRiccatiRecursionFactorizer::
~BackwardUnRiccatiRecursionFactorizer() {
}


inline void BackwardUnRiccatiRecursionFactorizer::factorizeKKTMatrix(
    const SplitRiccatiFactorization& riccati_next, const double dt, 
    SplitUnKKTMatrix& unkkt_matrix, SplitUnKKTResidual& unkkt_residual) {
  assert(dt > 0);
  // Factorize F
  unkkt_matrix.Qqq().noalias() += riccati_next.Pqq;
  unkkt_matrix.Qqv().noalias() += dt * riccati_next.Pqq;
  unkkt_matrix.Qqv().noalias() += riccati_next.Pqv;
  unkkt_matrix.Qvq() = unkkt_matrix.Qqv().transpose();
  unkkt_matrix.Qvv().noalias() += dt * dt * riccati_next.Pqq;
  unkkt_matrix.Qvv().noalias() += dt * riccati_next.Pqv;
  unkkt_matrix.Qvv().noalias() += dt * riccati_next.Pqv.transpose();
  unkkt_matrix.Qvv().noalias() += riccati_next.Pvv;
  // Factorize H
  unkkt_matrix.Qaq().transpose().noalias() += dt * riccati_next.Pqv;
  unkkt_matrix.Qav().transpose().noalias() += dt * dt * riccati_next.Pqv;
  unkkt_matrix.Qav().transpose().noalias() += dt * riccati_next.Pvv;
  // Factorize G
  unkkt_matrix.Qaa().noalias() += dt * dt * riccati_next.Pvv;
  // Factorize vector term
  unkkt_residual.la().noalias() += dt * riccati_next.Pqv.transpose() 
                                        * unkkt_residual.Fq();
  unkkt_residual.la().noalias() += dt * riccati_next.Pvv 
                                        * unkkt_residual.Fv();
  unkkt_residual.la().noalias() -= dt * riccati_next.sv;
}


inline void BackwardUnRiccatiRecursionFactorizer::factorizeRiccatiFactorization(
    const SplitRiccatiFactorization& riccati_next, 
    const SplitUnKKTMatrix& unkkt_matrix, 
    const SplitUnKKTResidual& unkkt_residual, 
    const LQRStateFeedbackPolicy& lqr_policy, const double dt, 
    SplitRiccatiFactorization& riccati) {
  assert(dt > 0);
  riccati.Pqq = unkkt_matrix.Qqq();
  riccati.Pqv = unkkt_matrix.Qqv();
  riccati.Pvv = unkkt_matrix.Qvv();
  GK_.noalias() = unkkt_matrix.Qaa() * lqr_policy.K; 
  riccati.Pqq.noalias() 
      -= lqr_policy.K.leftCols(dimv_).transpose() * GK_.leftCols(dimv_);
  riccati.Pqv.noalias() 
      -= lqr_policy.K.leftCols(dimv_).transpose() * GK_.rightCols(dimv_);
  riccati.Pvv.noalias() 
      -= lqr_policy.K.rightCols(dimv_).transpose() * GK_.rightCols(dimv_);
  riccati.Pvq = riccati.Pqv.transpose();
  // preserve the symmetry
  riccati.Pqq = 0.5 * (riccati.Pqq + riccati.Pqq.transpose()).eval();
  riccati.Pvv = 0.5 * (riccati.Pvv + riccati.Pvv.transpose()).eval();
  riccati.sq = riccati_next.sq;
  riccati.sq.noalias() -= riccati_next.Pqq * unkkt_residual.Fq();
  riccati.sq.noalias() -= riccati_next.Pqv * unkkt_residual.Fv();
  riccati.sv = riccati_next.sv;
  riccati.sv.noalias() += dt * riccati.sq;
  riccati.sv.noalias() -= riccati_next.Pqv.transpose() * unkkt_residual.Fq();
  riccati.sv.noalias() -= riccati_next.Pvv * unkkt_residual.Fv();
  riccati.sq.noalias() -= unkkt_residual.lq();
  riccati.sv.noalias() -= unkkt_residual.lv();
  riccati.sq.noalias() -= unkkt_matrix.Qaq().transpose() * lqr_policy.k;
  riccati.sv.noalias() -= unkkt_matrix.Qav().transpose() * lqr_policy.k;
}

} // namespace idocp

#endif // IDOCP_BACKWARD_UNRICCATI_RECURSION_FACTORIZER_HXX_ 