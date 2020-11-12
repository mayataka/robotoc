#ifndef IDOCP_RICCATI_FACTORIZER_HXX_
#define IDOCP_RICCATI_FACTORIZER_HXX_

#include "idocp/ocp/riccati_factorizer.hpp"

#include <cassert>

namespace idocp {

inline RiccatiFactorizer::RiccatiFactorizer(const Robot& robot) 
  : has_floating_base_(robot.has_floating_base()),
    has_state_constraint_(robot.max_point_contacts() > 0),
    dimv_(robot.dimv()),
    dimu_(robot.dimu()),
    llt_(robot.dimu()),
    K_(Eigen::MatrixXd::Zero(robot.dimu(), 2*robot.dimv())),
    k_(Eigen::VectorXd::Zero(robot.dimu())),
    AtPqq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    AtPqv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    AtPvq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    AtPvv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    BtPq_(Eigen::MatrixXd::Zero(robot.dimu(), robot.dimv())),
    BtPv_(Eigen::MatrixXd::Zero(robot.dimu(), robot.dimv())),
    GK_(Eigen::MatrixXd::Zero(robot.dimu(), 2*robot.dimv())),
    GinvBt_(Eigen::MatrixXd::Zero(robot.dimu(), robot.dimv())),
    BGinvBt_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    NApBKt_(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())) {
}


inline RiccatiFactorizer::RiccatiFactorizer() 
  : has_floating_base_(false),
    has_state_constraint_(false),
    dimv_(0),
    dimu_(0),
    llt_(),
    K_(),
    k_(),
    AtPqq_(),
    AtPqv_(),
    AtPvq_(),
    AtPvv_(),
    BtPq_(),
    BtPv_(),
    GK_(),
    GinvBt_(), 
    BGinvBt_(), 
    NApBKt_() {
}


inline RiccatiFactorizer::~RiccatiFactorizer() {
}


inline void RiccatiFactorizer::backwardRiccatiRecursion(
    const RiccatiFactorization& riccati_next, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual, const double dtau,
    RiccatiFactorization& riccati) {
  assert(dtau > 0);
  factorizeKKTMatrix(riccati_next, dtau, kkt_matrix, kkt_residual);
  computeFeedbackGainAndFeedforward(kkt_matrix, kkt_residual, riccati);
  factorizeRiccatiFactorization(riccati_next, kkt_matrix, kkt_residual, dtau, 
                                riccati);
}


inline void RiccatiFactorizer::forwardRiccatiRecursionParallel(
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual, 
    RiccatiFactorization& riccati) {
  kkt_matrix.Fxx().bottomRows(dimv_).noalias() += kkt_matrix.Fvu() * K_;
  GinvBt_.noalias() = llt_.solve(kkt_matrix.Fvu().transpose());
  BGinvBt_.noalias() = kkt_matrix.Fvu() * GinvBt_;
  kkt_residual.Fx().tail(dimv_).noalias() += kkt_matrix.Fvu() * k_;
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
    const RiccatiFactorization& riccati, SplitDirection& d, 
    const Eigen::MatrixBase<VectorType>& dx0) {
  d.dx().noalias() = riccati.Pi * dx0;
  d.dx().noalias() += riccati.pi;
  if (has_state_constraint_) {
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
  if (has_state_constraint_) {
    d.dlmdgmm().noalias() += riccati.n;
  }
}


inline void RiccatiFactorizer::computeControlInputDirection(
    const RiccatiFactorization& riccati, SplitDirection& d) const {
  d.du().noalias() = K_ * d.dx();
  d.du().noalias() += k_;
  if (has_state_constraint_) {
    d.du().noalias() -= GinvBt_ * riccati.n;
  }
}


template <typename MatrixType>
inline void RiccatiFactorizer::getStateFeedbackGain(
    const Eigen::MatrixBase<MatrixType>& K) {
  assert(K.rows() == dimu_);
  assert(K.cols() == 2*dimv_);
  const_cast<Eigen::MatrixBase<MatrixType>&> (K) = K_;
}


template <typename MatrixType1, typename MatrixType2>
inline void RiccatiFactorizer::getStateFeedbackGain(
    const Eigen::MatrixBase<MatrixType1>& Kq, 
    const Eigen::MatrixBase<MatrixType2>& Kv) {
  assert(Kq.rows() == dimu_);
  assert(Kq.cols() == dimv_);
  assert(Kv.rows() == dimu_);
  assert(Kv.cols() == dimv_);
  const_cast<Eigen::MatrixBase<MatrixType1>&> (Kq) = K_.leftCols(dimv_);
  const_cast<Eigen::MatrixBase<MatrixType2>&> (Kv) = K_.rightCols(dimv_);
}


inline void RiccatiFactorizer::factorizeKKTMatrix(
    const RiccatiFactorization& riccati_next, const double dtau, 
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  assert(dtau > 0);
  if (has_floating_base_) {
    AtPqq_.noalias() = kkt_matrix.Fqq().transpose() * riccati_next.Pqq;
    AtPqq_.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.Pqv.transpose();
    AtPqv_.noalias() = kkt_matrix.Fqq().transpose() * riccati_next.Pqv;
    AtPqv_.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.Pvv;
    AtPvq_ = dtau * riccati_next.Pqq;
    AtPvq_.noalias() += kkt_matrix.Fvv().transpose() * riccati_next.Pqv.transpose();
    AtPvv_ = dtau * riccati_next.Pqv;
    AtPvv_.noalias() += kkt_matrix.Fvv().transpose() * riccati_next.Pvv;
  }
  else {
    AtPqq_ = riccati_next.Pqq;
    AtPqq_.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.Pqv.transpose();
    AtPqv_ = riccati_next.Pqv;
    AtPqv_.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.Pvv;
    AtPvq_ = dtau * riccati_next.Pqq;
    AtPvq_.noalias() += kkt_matrix.Fvv().transpose() * riccati_next.Pqv.transpose();
    AtPvv_ = dtau * riccati_next.Pqv;
    AtPvv_.noalias() += kkt_matrix.Fvv().transpose() * riccati_next.Pvv;
  }
  BtPq_.noalias() = kkt_matrix.Fvu().transpose() * riccati_next.Pqv.transpose();
  BtPv_.noalias() = kkt_matrix.Fvu().transpose() * riccati_next.Pvv;
  // Factorize F
  if (has_floating_base_) {
    kkt_matrix.Qqq().noalias() += AtPqq_ * kkt_matrix.Fqq();
    kkt_matrix.Qqq().noalias() += AtPqv_ * kkt_matrix.Fvq();
    kkt_matrix.Qqv().noalias() += dtau * AtPqq_;
    kkt_matrix.Qqv().noalias() += AtPqv_ * kkt_matrix.Fvv();
    kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
    kkt_matrix.Qvv().noalias() += dtau * AtPvq_;
    kkt_matrix.Qvv().noalias() += AtPvv_ * kkt_matrix.Fvv();
  }
  else {
    kkt_matrix.Qqq().noalias() += AtPqq_;
    kkt_matrix.Qqq().noalias() += AtPqv_ * kkt_matrix.Fvq();
    kkt_matrix.Qqv().noalias() += dtau * AtPqq_;
    kkt_matrix.Qqv().noalias() += AtPqv_ * kkt_matrix.Fvv();
    kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
    kkt_matrix.Qvv().noalias() += dtau * AtPvq_;
    kkt_matrix.Qvv().noalias() += AtPvv_ * kkt_matrix.Fvv();
  }
  // Factorize H
  kkt_matrix.Qqu().noalias() += AtPqv_ * kkt_matrix.Fvu();
  kkt_matrix.Qvu().noalias() += AtPvv_ * kkt_matrix.Fvu();
  
  kkt_matrix.Quu().noalias() += BtPv_ * kkt_matrix.Fvu();
  // Factorize vector term
  kkt_residual.lu().noalias() += BtPq_ * kkt_residual.Fq();
  kkt_residual.lu().noalias() += BtPv_ * kkt_residual.Fv();
  kkt_residual.lu().noalias() -= kkt_matrix.Fvu().transpose() * riccati_next.sv;
}


inline void RiccatiFactorizer::computeFeedbackGainAndFeedforward(
    const KKTMatrix& kkt_matrix, const KKTResidual& kkt_residual, 
    RiccatiFactorization& riccati) {
  llt_.compute(kkt_matrix.Quu());
  assert(llt_.info() == Eigen::Success);
  K_ = - llt_.solve(kkt_matrix.Qxu().transpose());
  k_ = - llt_.solve(kkt_residual.lu());
}


inline void RiccatiFactorizer::factorizeRiccatiFactorization(
    const RiccatiFactorization& riccati_next, const KKTMatrix& kkt_matrix, 
    const KKTResidual& kkt_residual, const double dtau,
    RiccatiFactorization& riccati) {
  assert(dtau > 0);
  riccati.Pqq = kkt_matrix.Qqq();
  riccati.Pqv = kkt_matrix.Qqv();
  riccati.Pvv = kkt_matrix.Qvv();
  GK_.noalias() = kkt_matrix.Quu() * K_; 
  riccati.Pqq.noalias() -= riccati.Kq().transpose() * GK_.leftCols(dimv_);
  riccati.Pqv.noalias() -= riccati.Kq().transpose() * GK_.rightCols(dimv_);
  riccati.Pvv.noalias() -= riccati.Kv().transpose() * GK_.rightCols(dimv_);
  riccati.Pvq = riccati.Pqv.transpose();
  // preserve the symmetry
  riccati.Pqq = 0.5 * (riccati.Pqq + riccati.Pqq.transpose()).eval();
  riccati.Pvv = 0.5 * (riccati.Pvv + riccati.Pvv.transpose()).eval();
  if (has_floating_base_) {
    riccati.sq.noalias() = kkt_matrix.Fqq().transpose() * riccati_next.sq;
    riccati.sq.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.sv;
    riccati.sv.noalias() = dtau * riccati_next.sq;
    riccati.sv.noalias() += kkt_matrix.Fvv().transpose() * riccati_next.sv;
  }
  else {
    riccati.sq.noalias() = riccati_next.sq;
    riccati.sq.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.sv;
    riccati.sv.noalias() = dtau * riccati_next.sq;
    riccati.sv.noalias() += kkt_matrix.Fvv().transpose() * riccati_next.sv;
  }
  riccati.sq.noalias() -= AtPqq_ * kkt_residual.Fq();
  riccati.sq.noalias() -= AtPqv_ * kkt_residual.Fv();
  riccati.sv.noalias() -= AtPvq_ * kkt_residual.Fq();
  riccati.sv.noalias() -= AtPvv_ * kkt_residual.Fv();
  riccati.sq.noalias() -= kkt_residual.lq();
  riccati.sv.noalias() -= kkt_residual.lv();
  riccati.sq.noalias() -= kkt_matrix.Qqu() * k_;
  riccati.sv.noalias() -= kkt_matrix.Qvu() * k_;
}

} // namespace idocp

#endif // IDOCP_RICCATI_FACTORIZER_HXX_