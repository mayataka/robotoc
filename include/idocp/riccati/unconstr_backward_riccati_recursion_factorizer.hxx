#ifndef IDOCP_UNCONSTR_BACKWARD_RICCATI_RECURSION_FACTORIZER_HXX_
#define IDOCP_UNCONSTR_BACKWARD_RICCATI_RECURSION_FACTORIZER_HXX_

#include "idocp/riccati/unconstr_backward_riccati_recursion_factorizer.hpp"

#include <cassert>


namespace idocp {

inline UnconstrBackwardRiccatiRecursionFactorizer::
UnconstrBackwardRiccatiRecursionFactorizer(const Robot& robot) 
  : dimv_(robot.dimv()),
    GK_(Eigen::MatrixXd::Zero(robot.dimv(), 2*robot.dimv())),
    Pqq_tmp_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Pvv_tmp_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())) {
}


inline UnconstrBackwardRiccatiRecursionFactorizer::
UnconstrBackwardRiccatiRecursionFactorizer() 
  : dimv_(0),
    GK_(),
    Pqq_tmp_(),
    Pvv_tmp_() {
}


inline UnconstrBackwardRiccatiRecursionFactorizer::
~UnconstrBackwardRiccatiRecursionFactorizer() {
}


inline void UnconstrBackwardRiccatiRecursionFactorizer::factorizeKKTMatrix(
    const SplitRiccatiFactorization& riccati_next, const double dt, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  // Factorize F
  kkt_matrix.Qqq().noalias() += riccati_next.Pqq;
  kkt_matrix.Qqv().noalias() += dt * riccati_next.Pqq;
  kkt_matrix.Qqv().noalias() += riccati_next.Pqv;
  kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
  kkt_matrix.Qvv().noalias() += (dt*dt) * riccati_next.Pqq;
  kkt_matrix.Qvv().noalias() += dt * riccati_next.Pqv;
  kkt_matrix.Qvv().noalias() += dt * riccati_next.Pqv.transpose();
  kkt_matrix.Qvv().noalias() += riccati_next.Pvv;
  // Factorize H
  kkt_matrix.Qqu().noalias() += dt * riccati_next.Pqv;
  kkt_matrix.Qvu().noalias() += (dt*dt) * riccati_next.Pqv;
  kkt_matrix.Qvu().noalias() += dt * riccati_next.Pvv;
  // Factorize G
  kkt_matrix.Qaa.noalias()   += (dt*dt) * riccati_next.Pvv;
  // Factorize vector term
  kkt_residual.la.noalias() += dt * riccati_next.Pqv.transpose() * kkt_residual.Fq();
  kkt_residual.la.noalias() += dt * riccati_next.Pvv * kkt_residual.Fv();
  kkt_residual.la.noalias() -= dt * riccati_next.sv;
}


inline void UnconstrBackwardRiccatiRecursionFactorizer::
factorizeRiccatiFactorization(const SplitRiccatiFactorization& riccati_next, 
                              const SplitKKTMatrix& kkt_matrix, 
                              const SplitKKTResidual& kkt_residual, 
                              const LQRPolicy& lqr_policy, const double dt, 
                              SplitRiccatiFactorization& riccati) {
  assert(dt > 0);
  riccati.Pqq = kkt_matrix.Qqq();
  riccati.Pqv = kkt_matrix.Qqv();
  riccati.Pvv = kkt_matrix.Qvv();
  GK_.noalias() = kkt_matrix.Qaa * lqr_policy.K; 
  riccati.Pqq.noalias() 
      -= lqr_policy.K.leftCols(dimv_).transpose() * GK_.leftCols(dimv_);
  riccati.Pqv.noalias() 
      -= lqr_policy.K.leftCols(dimv_).transpose() * GK_.rightCols(dimv_);
  riccati.Pvv.noalias() 
      -= lqr_policy.K.rightCols(dimv_).transpose() * GK_.rightCols(dimv_);
  riccati.Pvq = riccati.Pqv.transpose();
  // preserve the symmetry
  Pqq_tmp_ = riccati.Pqq + riccati.Pqq.transpose();
  Pvv_tmp_ = riccati.Pvv + riccati.Pvv.transpose();
  riccati.Pqq = 0.5 * Pqq_tmp_;
  riccati.Pvv = 0.5 * Pvv_tmp_;
  riccati.sq = riccati_next.sq;
  riccati.sq.noalias() -= riccati_next.Pqq * kkt_residual.Fq();
  riccati.sq.noalias() -= riccati_next.Pqv * kkt_residual.Fv();
  riccati.sv = riccati_next.sv;
  riccati.sv.noalias() += dt * riccati.sq;
  riccati.sv.noalias() -= riccati_next.Pqv.transpose() * kkt_residual.Fq();
  riccati.sv.noalias() -= riccati_next.Pvv * kkt_residual.Fv();
  riccati.sq.noalias() -= kkt_residual.lq();
  riccati.sv.noalias() -= kkt_residual.lv();
  riccati.sq.noalias() -= kkt_matrix.Qqu() * lqr_policy.k;
  riccati.sv.noalias() -= kkt_matrix.Qvu() * lqr_policy.k;
}

} // namespace idocp

#endif // IDOCP_UNCONSTR_BACKWARD_RICCATI_RECURSION_FACTORIZER_HXX_ 