#ifndef IDOCP_SPLIT_UNRICCATI_FACTORIZER_HXX_ 
#define IDOCP_SPLIT_UNRICCATI_FACTORIZER_HXX_

#include "idocp/unocp/split_unriccati_factorizer.hpp"

#include <cassert>

namespace idocp {

inline SplitUnRiccatiFactorizer::SplitUnRiccatiFactorizer(const Robot& robot) 
  : dimv_(robot.dimv()),
    llt_(robot.dimv()),
    lqr_policy_(robot),
    backward_recursion_(robot) {
}


inline SplitUnRiccatiFactorizer::SplitUnRiccatiFactorizer() 
  : dimv_(0),
    llt_(),
    lqr_policy_(),
    backward_recursion_() {
}


inline SplitUnRiccatiFactorizer::~SplitUnRiccatiFactorizer() {
}


inline void SplitUnRiccatiFactorizer::backwardRiccatiRecursion(
    const SplitRiccatiFactorization& riccati_next, const double dtau, 
    SplitUnKKTMatrix& unkkt_matrix, SplitUnKKTResidual& unkkt_residual,  
    SplitRiccatiFactorization& riccati) {
  assert(dtau >= 0);
  backward_recursion_.factorizeKKTMatrix(riccati_next, dtau, unkkt_matrix, 
                                         unkkt_residual);
  llt_.compute(unkkt_matrix.Qaa());
  assert(llt_.info() == Eigen::Success);
  lqr_policy_.K = - llt_.solve(unkkt_matrix.Qax());
  lqr_policy_.k = - llt_.solve(unkkt_residual.la());
  assert(!lqr_policy_.K.hasNaN());
  assert(!lqr_policy_.k.hasNaN());
  backward_recursion_.factorizeRiccatiFactorization(riccati_next, unkkt_matrix, 
                                                    unkkt_residual, lqr_policy_,
                                                    dtau, riccati);
}


inline void SplitUnRiccatiFactorizer::forwardRiccatiRecursion(
    const SplitUnKKTResidual& unkkt_residual, SplitDirection& d, 
    const double dtau, SplitDirection& d_next) const {
  assert(dtau >= 0);
  d.da().noalias() = lqr_policy_.K * d.dx() + lqr_policy_.k;
  d_next.dx() = unkkt_residual.Fx() + d.dx();
  d_next.dq().noalias() += dtau * d.dv();
  d_next.dv().noalias() += dtau * d.da();
}


inline void SplitUnRiccatiFactorizer::computeCostateDirection(
    const SplitRiccatiFactorization& riccati, SplitDirection& d) {
  d.dlmd().noalias()  = riccati.Pqq * d.dq();
  d.dlmd().noalias() += riccati.Pqv * d.dv();
  d.dlmd().noalias() -= riccati.sq;
  d.dgmm().noalias()  = riccati.Pqv.transpose() * d.dq();
  d.dgmm().noalias() += riccati.Pvv * d.dv();
  d.dgmm().noalias() -= riccati.sv;
}


template <typename MatrixType>
inline void SplitUnRiccatiFactorizer::getStateFeedbackGain(
    const Eigen::MatrixBase<MatrixType>& K) const {
  assert(K.rows() == dimv_);
  assert(K.cols() == 2*dimv_);
  const_cast<Eigen::MatrixBase<MatrixType>&> (K) = lqr_policy_.K;
}


template <typename MatrixType1, typename MatrixType2>
inline void SplitUnRiccatiFactorizer::getStateFeedbackGain(
    const Eigen::MatrixBase<MatrixType1>& Kq, 
    const Eigen::MatrixBase<MatrixType2>& Kv) const {
  assert(Kq.rows() == dimv_);
  assert(Kq.cols() == dimv_);
  assert(Kv.rows() == dimv_);
  assert(Kv.cols() == dimv_);
  const_cast<Eigen::MatrixBase<MatrixType1>&> (Kq) 
      = lqr_policy_.K.leftCols(dimv_);
  const_cast<Eigen::MatrixBase<MatrixType2>&> (Kv) 
      = lqr_policy_.K.rightCols(dimv_);
}

} // namespace idocp

#endif // IDOCP_SPLIT_UNRICCATI_FACTORIZER_HXX_ 