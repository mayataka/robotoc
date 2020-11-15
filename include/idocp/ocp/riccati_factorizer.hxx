#ifndef IDOCP_RICCATI_FACTORIZER_HXX_
#define IDOCP_RICCATI_FACTORIZER_HXX_

#include "idocp/ocp/riccati_factorizer.hpp"

#include <cassert>

namespace idocp {

inline RiccatiFactorizer::RiccatiFactorizer(const Robot& robot) 
  : has_floating_base_(robot.has_floating_base()),
    dimv_(robot.dimv()),
    dimu_(robot.dimu()),
    llt_(robot.dimu()),
    lqr_policy_(robot),
    backward_recursion_(robot),
    forward_recursion_(robot),
    GinvBt_(Eigen::MatrixXd::Zero(robot.dimu(), robot.dimv())),
    BGinvBt_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())) {
}


inline RiccatiFactorizer::RiccatiFactorizer() 
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


inline RiccatiFactorizer::~RiccatiFactorizer() {
}


inline void RiccatiFactorizer::backwardRiccatiRecursion(
    const RiccatiFactorization& riccati_next, const double dtau, 
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual, 
    RiccatiFactorization& riccati) {
  assert(dtau > 0);
  backward_recursion_.factorizeKKTMatrix(riccati_next, dtau, kkt_matrix, 
                                         kkt_residual);
  llt_.compute(kkt_matrix.Quu());
  assert(llt_.info() == Eigen::Success);
  lqr_policy_.K = - llt_.solve(kkt_matrix.Qxu().transpose());
  lqr_policy_.k = - llt_.solve(kkt_residual.lu());
  backward_recursion_.factorizeRiccatiFactorization(riccati_next, kkt_matrix, 
                                                    kkt_residual, lqr_policy_,
                                                    dtau, riccati);
}


inline void RiccatiFactorizer::forwardRiccatiRecursionParallel(
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual, 
    const bool exist_state_constraint) {
  kkt_matrix.Fxx().bottomRows(dimv_).noalias() 
      += kkt_matrix.Fvu() * lqr_policy_.K;
  kkt_residual.Fx().tail(dimv_).noalias() += kkt_matrix.Fvu() * lqr_policy_.k;
  if (exist_state_constraint) {
    GinvBt_ = llt_.solve(kkt_matrix.Fvu().transpose());
    BGinvBt_.noalias() = kkt_matrix.Fvu() * GinvBt_;
  }
}


inline void RiccatiFactorizer::forwardRiccatiRecursionSerialInitial(
    const RiccatiFactorization& riccati) {
  assert(riccati.Pi.isIdentity()); // assume riccati.Pi is already a identity matrix
  assert(riccati.pi.isZero()); // assume riccati.pi is already a zero vector  
  assert(riccati.N.isZero()); // assume riccati.N is already a zero matrix.
}


inline void RiccatiFactorizer::forwardRiccatiRecursionSerial(
    const RiccatiFactorization& riccati, const KKTMatrix& kkt_matrix, 
    const KKTResidual& kkt_residual, const double dtau,
    RiccatiFactorization& riccati_next, 
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
inline void RiccatiFactorizer::backwardStateConstraintFactorization(
    const KKTMatrix& kkt_matrix, const Eigen::MatrixBase<MatrixType1>& T_next,
    const double dtau, const Eigen::MatrixBase<MatrixType2>& T) const {
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


inline void RiccatiFactorizer::forwardRiccatiRecursion(
    const KKTMatrix& kkt_matrix, const KKTResidual& kkt_residual, 
    SplitDirection& d, const double dtau, 
    SplitDirection& d_next) const {
  assert(dtau > 0);
  d.du().noalias() = lqr_policy_.K * d.dx();
  d.du().noalias() += lqr_policy_.k;
  if (has_floating_base_) {
    d_next.dq().noalias() = kkt_matrix.Fqq() * d.dq();
    d_next.dq().noalias() += dtau * d.dv() + kkt_residual.Fq();
  }
  else {
    d_next.dq().noalias() = d.dq() + dtau * d.dv() + kkt_residual.Fq();
  }
  d_next.dv().noalias() = kkt_matrix.Fvq() * d.dq();
  d_next.dv().noalias() += kkt_matrix.Fvv() * d.dv();
  // d_next.dv().noalias() += kkt_matrix.Fvu() * d.du();
  d_next.dv().noalias() += kkt_residual.Fv();
}


template <typename VectorType>
inline void RiccatiFactorizer::computeStateDirection(
    const RiccatiFactorization& riccati, 
    const Eigen::MatrixBase<VectorType>& dx0, SplitDirection& d,
    const bool exist_state_constraint) {
  d.dx().noalias() = riccati.Pi * dx0;
  d.dx().noalias() += riccati.pi;
  if (exist_state_constraint) {
    d.dx().noalias() -= riccati.N * riccati.n;
  }
}


inline void RiccatiFactorizer::computeCostateDirection(
    const RiccatiFactorization& riccati, SplitDirection& d,
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


inline void RiccatiFactorizer::computeControlInputDirection(
    const RiccatiFactorization& riccati, SplitDirection& d,
    const bool exist_state_constraint) const {
  d.du().noalias() = lqr_policy_.K * d.dx();
  d.du().noalias() += lqr_policy_.k;
  if (exist_state_constraint) {
    d.du().noalias() -= GinvBt_ * riccati.n.tail(dimv_);
  }
}


template <typename MatrixType>
inline void RiccatiFactorizer::getStateFeedbackGain(
    const Eigen::MatrixBase<MatrixType>& K) {
  assert(K.rows() == dimu_);
  assert(K.cols() == 2*dimv_);
  const_cast<Eigen::MatrixBase<MatrixType>&> (K) = lqr_policy_.K;
}


template <typename MatrixType1, typename MatrixType2>
inline void RiccatiFactorizer::getStateFeedbackGain(
    const Eigen::MatrixBase<MatrixType1>& Kq, 
    const Eigen::MatrixBase<MatrixType2>& Kv) {
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

#endif // IDOCP_RICCATI_FACTORIZER_HXX_