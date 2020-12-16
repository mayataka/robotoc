#ifndef IDOCP_SPLIT_RICCATI_FACTORIZER_HXX_ 
#define IDOCP_SPLIT_RICCATI_FACTORIZER_HXX_

#include "idocp/ocp/split_riccati_factorizer.hpp"

#include <cassert>

namespace idocp {

inline SplitRiccatiFactorizer::SplitRiccatiFactorizer(const Robot& robot) 
  : has_floating_base_(robot.hasFloatingBase()),
    dimv_(robot.dimv()),
    dimu_(robot.dimu()),
    llt_(robot.dimu()),
    lqr_policy_(robot),
    backward_recursion_(robot),
    forward_recursion_(robot),
    GinvBt_(Eigen::MatrixXd::Zero(robot.dimu(), robot.dimv())),
    BGinvBt_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())) {
}


inline SplitRiccatiFactorizer::SplitRiccatiFactorizer() 
  : has_floating_base_(false),
    dimv_(0),
    dimu_(0),
    llt_(),
    lqr_policy_(),
    backward_recursion_(),
    forward_recursion_(),
    GinvBt_(), 
    BGinvBt_() {
}


inline SplitRiccatiFactorizer::~SplitRiccatiFactorizer() {
}


inline void SplitRiccatiFactorizer::backwardRiccatiRecursion(
    const SplitRiccatiFactorization& riccati_next, const double dtau, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual, 
    SplitRiccatiFactorization& riccati) {
  assert(dtau > 0);
  backward_recursion_.factorizeKKTMatrix(riccati_next, dtau, kkt_matrix, 
                                         kkt_residual);
  llt_.compute(kkt_matrix.Quu());
  assert(llt_.info() == Eigen::Success);
  lqr_policy_.K = - llt_.solve(kkt_matrix.Qxu().transpose());
  lqr_policy_.k = - llt_.solve(kkt_residual.lu());
  assert(!lqr_policy_.K.hasNaN());
  assert(!lqr_policy_.k.hasNaN());
  backward_recursion_.factorizeRiccatiFactorization(riccati_next, kkt_matrix, 
                                                    kkt_residual, lqr_policy_,
                                                    dtau, riccati);
}


inline void SplitRiccatiFactorizer::forwardRiccatiRecursionParallel(
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual, 
    const bool exist_state_constraint) {
  kkt_matrix.Fxx().bottomRows(dimv_).noalias() 
      += kkt_matrix.Fvu() * lqr_policy_.K;
  kkt_residual.Fx().tail(dimv_).noalias() += kkt_matrix.Fvu() * lqr_policy_.k;
  if (exist_state_constraint) {
    GinvBt_ = llt_.solve(kkt_matrix.Fvu().transpose());
    BGinvBt_.noalias() = kkt_matrix.Fvu() * GinvBt_;
  }
}


inline void SplitRiccatiFactorizer::forwardStateConstraintFactorizationInitial(
    const SplitRiccatiFactorization& riccati) {
  assert(riccati.Pi.isIdentity()); // Checks riccati.Pi is a identity matrix.
  assert(riccati.pi.isZero()); // Checks riccati.pi is a zero vector.
  assert(riccati.N.isZero()); // Checks riccati.N is a zero matrix.
}


inline void SplitRiccatiFactorizer::forwardStateConstraintFactorization(
    const SplitRiccatiFactorization& riccati, 
    const SplitKKTMatrix& kkt_matrix, const SplitKKTResidual& kkt_residual, 
    const double dtau, SplitRiccatiFactorization& riccati_next, 
    const bool exist_state_constraint) {
  assert(dtau > 0);
  forward_recursion_.factorizeStateTransition(riccati, kkt_matrix, kkt_residual, 
                                              dtau, riccati_next);
  if (exist_state_constraint) {
    forward_recursion_.factorizeStateConstraintFactorization(riccati, 
                                                             kkt_matrix, dtau, 
                                                             riccati_next);
    riccati_next.N.bottomRightCorner(dimv_, dimv_).noalias() += BGinvBt_;
  }
}


template <typename MatrixType1, typename MatrixType2>
inline void SplitRiccatiFactorizer::backwardStateConstraintFactorization(
    const Eigen::MatrixBase<MatrixType1>& T_next, 
    const SplitKKTMatrix& kkt_matrix, const double dtau, 
    const Eigen::MatrixBase<MatrixType2>& T) const {
  assert(T_next.rows() == T.rows());
  assert(T_next.rows() == T.rows());
  assert(dtau > 0);
  if (has_floating_base_) {
    const_cast<Eigen::MatrixBase<MatrixType2>&> (T).topRows(dimv_).noalias() 
        = kkt_matrix.Fqq().transpose() * T_next.topRows(dimv_);
  }
  else {
    const_cast<Eigen::MatrixBase<MatrixType2>&> (T).topRows(dimv_) 
        = T_next.topRows(dimv_);
  }
  const_cast<Eigen::MatrixBase<MatrixType2>&> (T).topRows(dimv_).noalias()
      += kkt_matrix.Fvq().transpose() * T_next.bottomRows(dimv_);
  const_cast<Eigen::MatrixBase<MatrixType2>&> (T).bottomRows(dimv_)
      = dtau * T_next.topRows(dimv_);
  const_cast<Eigen::MatrixBase<MatrixType2>&> (T).bottomRows(dimv_).noalias() 
      +=  kkt_matrix.Fvv().transpose() * T_next.bottomRows(dimv_);
}


template <typename SplitDirectionType>
inline void SplitRiccatiFactorizer::forwardRiccatiRecursion(
    const SplitKKTMatrix& kkt_matrix, const SplitKKTResidual& kkt_residual, 
    const SplitRiccatiFactorization& riccati_next, const SplitDirection& d, 
    const double dtau, SplitDirectionType& d_next, 
    const bool exist_state_constraint) const {
  assert(dtau > 0);
  d_next.dx() = kkt_residual.Fx();
  if (has_floating_base_) {
    d_next.dq().noalias() += kkt_matrix.Fqq() * d.dq();
    d_next.dq().noalias() += dtau * d.dv();
  }
  else {
    d_next.dq().noalias() += d.dq() + dtau * d.dv();
  }
  d_next.dv().noalias() += kkt_matrix.Fvq() * d.dq();
  d_next.dv().noalias() += kkt_matrix.Fvv() * d.dv();
  if (exist_state_constraint) {
    d_next.dv().noalias() -= BGinvBt_ * riccati_next.n.tail(dimv_);
  }
}


inline void SplitRiccatiFactorizer::computeCostateDirection(
    const SplitRiccatiFactorization& riccati, SplitDirection& d,
    const bool exist_state_constraint) {
  d.dlmd().noalias() = riccati.Pqq * d.dq();
  d.dlmd().noalias() += riccati.Pqv * d.dv();
  d.dlmd().noalias() -= riccati.sq;
  d.dgmm().noalias() = riccati.Pqv.transpose() * d.dq();
  d.dgmm().noalias() += riccati.Pvv * d.dv();
  d.dgmm().noalias() -= riccati.sv;
  if (exist_state_constraint) {
    d.dlmdgmm().noalias() += riccati.n;
  }
}


inline void SplitRiccatiFactorizer::computeControlInputDirection(
    const SplitRiccatiFactorization& riccati_next, SplitDirection& d,
    const bool exist_state_constraint) const {
  d.du().noalias() = lqr_policy_.K * d.dx();
  d.du().noalias() += lqr_policy_.k;
  if (exist_state_constraint) {
    d.du().noalias() -= GinvBt_ * riccati_next.n.tail(dimv_);
  }
}


template <typename MatrixType>
inline void SplitRiccatiFactorizer::getStateFeedbackGain(
    const Eigen::MatrixBase<MatrixType>& K) const {
  assert(K.rows() == dimu_);
  assert(K.cols() == 2*dimv_);
  const_cast<Eigen::MatrixBase<MatrixType>&> (K) = lqr_policy_.K;
}


template <typename MatrixType1, typename MatrixType2>
inline void SplitRiccatiFactorizer::getStateFeedbackGain(
    const Eigen::MatrixBase<MatrixType1>& Kq, 
    const Eigen::MatrixBase<MatrixType2>& Kv) const {
  assert(Kq.rows() == dimu_);
  assert(Kq.cols() == dimv_);
  assert(Kv.rows() == dimu_);
  assert(Kv.cols() == dimv_);
  const_cast<Eigen::MatrixBase<MatrixType1>&> (Kq) 
      = lqr_policy_.K.leftCols(dimv_);
  const_cast<Eigen::MatrixBase<MatrixType2>&> (Kv) 
      = lqr_policy_.K.rightCols(dimv_);
}

} // namespace idocp

#endif // IDOCP_SPLIT_RICCATI_FACTORIZER_HXX_ 