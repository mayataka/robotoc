#ifndef IDOCP_RICCATI_FACTORIZER_HXX_
#define IDOCP_RICCATI_FACTORIZER_HXX_

#include "idocp/ocp/riccati_factorizer.hpp"

#include <cassert>

namespace idocp {

inline RiccatiFactorizer::RiccatiFactorizer(const Robot& robot) 
  : has_floating_base_(robot.has_floating_base()),
    exist_state_constraint_(robot.max_point_contacts() > 0),
    dimv_(robot.dimv()),
    dimu_(robot.dimu()),
    llt_(robot.dimu()),
    lqr_policy_(robot),
    backward_recursion_(robot),
    GinvBt_(Eigen::MatrixXd::Zero(robot.dimu(), robot.dimv())),
    BGinvBt_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    NApBKt_(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())) {
}


inline RiccatiFactorizer::RiccatiFactorizer() 
  : has_floating_base_(false),
    exist_state_constraint_(false),
    dimv_(0),
    dimu_(0),
    llt_(),
    lqr_policy_(),
    backward_recursion_(),
    GinvBt_(), 
    BGinvBt_(), 
    NApBKt_() {
}


inline RiccatiFactorizer::~RiccatiFactorizer() {
}


inline void RiccatiFactorizer::setExistStateConstraint(
    const bool exist_state_constraint) {
  exist_state_constraint_ = exist_state_constraint;
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
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  kkt_matrix.Fxx().bottomRows(dimv_).noalias() += kkt_matrix.Fvu() * lqr_policy_.K;
  GinvBt_.noalias() = llt_.solve(kkt_matrix.Fvu().transpose());
  BGinvBt_.noalias() = kkt_matrix.Fvu() * GinvBt_;
  kkt_residual.Fx().tail(dimv_).noalias() += kkt_matrix.Fvu() * lqr_policy_.k;
}


inline void RiccatiFactorizer::forwardRiccatiRecursionSerialInitial(
    const KKTMatrix& kkt_matrix, const KKTResidual& kkt_residual, 
    RiccatiFactorization& riccati) {
  assert(riccati.Pi.isIdentity()); // assume riccati.Pi is already a identity matrix
  assert(riccati.pi.isZero()); // assume riccati.pi is already a zero vector  
  assert(riccati.N.isZero()); // assume riccati.N is already a zero matrix.
}


inline void RiccatiFactorizer::forwardRiccatiRecursionSerial(
    const RiccatiFactorization& riccati, const KKTMatrix& kkt_matrix, 
    const KKTResidual& kkt_residual, RiccatiFactorization& riccati_next) {
  riccati_next.Pi.noalias() = kkt_matrix.Fxx() * riccati.Pi;
  riccati_next.pi = kkt_residual.Fx();
  riccati_next.pi.noalias() += kkt_matrix.Fxx() * riccati.pi;
  NApBKt_.noalias() = riccati.N * kkt_matrix.Fxx().transpose();
  riccati_next.N.noalias() = kkt_matrix.Fxx() * NApBKt_;
  riccati_next.N.bottomRightCorner(dimv_, dimv_).noalias() += BGinvBt_;
}


template <typename SplitDirectionType>
inline void RiccatiFactorizer::forwardRiccatiRecursion(
    const KKTMatrix& kkt_matrix, const KKTResidual& kkt_residual, 
    const SplitDirection& d, const double dtau, 
    SplitDirectionType& d_next) const {
  assert(dtau > 0);
  if (has_floating_base_) {
    d_next.dq().noalias() = kkt_matrix.Fqq() * d.dq();
    d_next.dq().noalias() += dtau * d.dv() + kkt_residual.Fq();
  }
  else {
    d_next.dq().noalias() = d.dq() + dtau * d.dv() + kkt_residual.Fq();
  }
  d_next.dv().noalias() = kkt_matrix.Fvq() * d.dq();
  d_next.dv().noalias() += kkt_matrix.Fvv() * d.dv();
  d_next.dv().noalias() += kkt_matrix.Fvu() * d.du();
  d_next.dv().noalias() += kkt_residual.Fv();
}


template <typename VectorType>
inline void RiccatiFactorizer::computeStateDirection(
    const RiccatiFactorization& riccati, 
    const Eigen::MatrixBase<VectorType>& dx0, SplitDirection& d) {
  d.dx().noalias() = riccati.Pi * dx0;
  d.dx().noalias() += riccati.pi;
  if (exist_state_constraint_) {
    d.dx().noalias() -= riccati.N * riccati.n;
  }
}


inline void RiccatiFactorizer::computeCostateDirection(
    const RiccatiFactorization& riccati, SplitDirection& d) const {
  d.dlmd().noalias() = riccati.Pqq * d.dq();
  d.dlmd().noalias() += riccati.Pqv * d.dv();
  d.dlmd().noalias() -= riccati.sq;
  d.dgmm().noalias() = riccati.Pqv.transpose() * d.dq();
  d.dgmm().noalias() += riccati.Pvv * d.dv();
  d.dgmm().noalias() -= riccati.sv;
  if (exist_state_constraint_) {
    d.dlmdgmm().noalias() += riccati.n;
  }
}


inline void RiccatiFactorizer::computeControlInputDirection(
    const RiccatiFactorization& riccati, SplitDirection& d) const {
  d.du().noalias() = lqr_policy_.K * d.dx();
  d.du().noalias() += lqr_policy_.k;
  if (exist_state_constraint_) {
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