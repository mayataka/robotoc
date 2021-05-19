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
    llt_s_(),
    lqr_policy_(robot),
    backward_recursion_(robot) {
}


inline SplitRiccatiFactorizer::SplitRiccatiFactorizer() 
  : has_floating_base_(false),
    dimv_(0),
    dimu_(0),
    llt_(),
    llt_s_(),
    lqr_policy_(),
    backward_recursion_() {
}


inline SplitRiccatiFactorizer::~SplitRiccatiFactorizer() {
}


inline void SplitRiccatiFactorizer::backwardRiccatiRecursion(
    const SplitRiccatiFactorization& riccati_next, const double dt, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual, 
    SplitRiccatiFactorization& riccati) {
  assert(dt > 0);
  backward_recursion_.factorizeKKTMatrix(riccati_next, dt, kkt_matrix, 
                                         kkt_residual);
  llt_.compute(kkt_matrix.Quu);
  assert(llt_.info() == Eigen::Success);
  lqr_policy_.K.noalias() = - llt_.solve(kkt_matrix.Qxu.transpose());
  lqr_policy_.k.noalias() = - llt_.solve(kkt_residual.lu);
  assert(!lqr_policy_.K.hasNaN());
  assert(!lqr_policy_.k.hasNaN());
  backward_recursion_.factorizeRiccatiFactorization(riccati_next, kkt_matrix, 
                                                    kkt_residual, lqr_policy_,
                                                    dt, riccati);
}


inline void SplitRiccatiFactorizer::backwardRiccatiRecursion(
    const SplitRiccatiFactorization& riccati_next, const double dt, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual, 
    const SplitStateConstraintJacobian& jac, SplitRiccatiFactorization& riccati,
    SplitConstrainedRiccatiFactorization& c_riccati) {
  assert(dt > 0);
  backward_recursion_.factorizeKKTMatrix(riccati_next, dt, kkt_matrix, 
                                         kkt_residual);
  // Schur complement
  llt_.compute(kkt_matrix.Quu);
  assert(llt_.info() == Eigen::Success);
  c_riccati.setImpulseStatus(jac.dimi());
  c_riccati.Ginv.noalias() = llt_.solve(Eigen::MatrixXd::Identity(dimu_, dimu_));
  c_riccati.DGinv().transpose().noalias() = llt_.solve(jac.Phiu().transpose());
  c_riccati.S().noalias() = c_riccati.DGinv() * jac.Phiu().transpose();
  llt_s_.compute(c_riccati.S());
  assert(llt_s_.info() == Eigen::Success);
  c_riccati.SinvDGinv().noalias() = llt_s_.solve(c_riccati.DGinv());
  c_riccati.Ginv.noalias() -= c_riccati.SinvDGinv().transpose() * c_riccati.DGinv();
  lqr_policy_.K.noalias()  = - c_riccati.Ginv * kkt_matrix.Qxu.transpose();
  lqr_policy_.K.noalias() -= c_riccati.SinvDGinv().transpose() * jac.Phix();
  lqr_policy_.k.noalias()  = - c_riccati.Ginv * kkt_residual.lu;
  lqr_policy_.k.noalias() -= c_riccati.SinvDGinv().transpose() * kkt_residual.P();
  c_riccati.M().noalias()  = llt_s_.solve(jac.Phix());
  c_riccati.M().noalias() -= c_riccati.SinvDGinv() * kkt_matrix.Qxu.transpose();
  c_riccati.m().noalias()  = llt_s_.solve(kkt_residual.P());
  c_riccati.m().noalias() -= c_riccati.SinvDGinv() * kkt_residual.lu;
  assert(!lqr_policy_.K.hasNaN());
  assert(!lqr_policy_.k.hasNaN());
  assert(!c_riccati.M().hasNaN());
  assert(!c_riccati.m().hasNaN());
  backward_recursion_.factorizeRiccatiFactorization(riccati_next, kkt_matrix, 
                                                    kkt_residual, lqr_policy_,
                                                    dt, riccati);
  c_riccati.DtM.noalias()   = jac.Phiu().transpose() * c_riccati.M();
  c_riccati.KtDtM.noalias() = lqr_policy_.K.transpose() * c_riccati.DtM;
  riccati.Pqq.noalias() -= c_riccati.KtDtM.topLeftCorner(dimv_, dimv_);
  riccati.Pqq.noalias() -= c_riccati.KtDtM.topLeftCorner(dimv_, dimv_).transpose();
  riccati.Pqv.noalias() -= c_riccati.KtDtM.topRightCorner(dimv_, dimv_);
  riccati.Pqv.noalias() -= c_riccati.KtDtM.bottomLeftCorner(dimv_, dimv_).transpose();
  riccati.Pvq = riccati.Pqv.transpose();
  riccati.Pvv.noalias() -= c_riccati.KtDtM.bottomRightCorner(dimv_, dimv_);
  riccati.Pvv.noalias() -= c_riccati.KtDtM.bottomRightCorner(dimv_, dimv_).transpose();
  riccati.sq.noalias() -= jac.Phix().transpose().topRows(dimv_) * c_riccati.m();
  riccati.sv.noalias() -= jac.Phix().transpose().bottomRows(dimv_) * c_riccati.m();
}


template <typename SplitDirectionType>
inline void SplitRiccatiFactorizer::forwardRiccatiRecursion(
    const SplitKKTMatrix& kkt_matrix, const SplitKKTResidual& kkt_residual, 
    const double dt, SplitDirection& d, SplitDirectionType& d_next) const {
  assert(dt > 0);
  d.du.noalias()  = lqr_policy_.K * d.dx;
  d.du.noalias() += lqr_policy_.k;
  d_next.dx = kkt_residual.Fx;
  if (has_floating_base_) {
    d_next.dq().template head<6>().noalias() 
        += kkt_matrix.Fqq().template topLeftCorner<6, 6>() 
            * d.dq().template head<6>();
    d_next.dq().tail(dimv_-6).noalias() += d.dq().tail(dimv_-6);
    d_next.dq().template head<6>().noalias() 
        += kkt_matrix.Fqv().template topLeftCorner<6, 6>() 
            * d.dv().template head<6>();
    d_next.dq().tail(dimv_-6).noalias() += dt * d.dv().tail(dimv_-6);
  }
  else {
    d_next.dq().noalias() += d.dq();
    d_next.dq().noalias() += dt * d.dv();
  }
  d_next.dv().noalias() += kkt_matrix.Fvq() * d.dq();
  d_next.dv().noalias() += kkt_matrix.Fvv() * d.dv();
  d_next.dv().noalias() += kkt_matrix.Fvu * d.du;
}


inline void SplitRiccatiFactorizer::computeCostateDirection(
    const SplitRiccatiFactorization& riccati, SplitDirection& d) {
  d.dlmd().noalias()  = riccati.Pqq * d.dq();
  d.dlmd().noalias() += riccati.Pqv * d.dv();
  d.dlmd().noalias() -= riccati.sq;
  d.dgmm().noalias()  = riccati.Pqv.transpose() * d.dq();
  d.dgmm().noalias() += riccati.Pvv * d.dv();
  d.dgmm().noalias() -= riccati.sv;
}


inline void SplitRiccatiFactorizer::
computeLagrangeMultiplierDirection(
    const SplitConstrainedRiccatiFactorization& c_riccati, 
    SplitDirection& d) {
  d.dxi().noalias()  = c_riccati.M() * d.dx;
  d.dxi().noalias() += c_riccati.m();
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