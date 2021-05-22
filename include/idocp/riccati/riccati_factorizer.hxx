#ifndef IDOCP_RICCATI_FACTORIZER_HXX_ 
#define IDOCP_RICCATI_FACTORIZER_HXX_

#include "idocp/riccati/riccati_factorizer.hpp"

#include <cassert>

namespace idocp {

inline RiccatiFactorizer::RiccatiFactorizer(const Robot& robot) 
  : has_floating_base_(robot.hasFloatingBase()),
    dimv_(robot.dimv()),
    dimu_(robot.dimu()),
    llt_(robot.dimu()),
    llt_s_(),
    backward_recursion_(robot) {
}


inline RiccatiFactorizer::RiccatiFactorizer() 
  : has_floating_base_(false),
    dimv_(0),
    dimu_(0),
    llt_(),
    llt_s_(),
    backward_recursion_() {
}


inline RiccatiFactorizer::~RiccatiFactorizer() {
}


inline void RiccatiFactorizer::backwardRiccatiRecursion(
    const SplitRiccatiFactorization& riccati_next,  
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual, 
    SplitRiccatiFactorization& riccati, LQRPolicy& lqr_policy) {
  backward_recursion_.factorizeKKTMatrix(riccati_next, kkt_matrix, 
                                         kkt_residual);
  llt_.compute(kkt_matrix.Quu);
  assert(llt_.info() == Eigen::Success);
  lqr_policy.K.noalias() = - llt_.solve(kkt_matrix.Qxu.transpose());
  lqr_policy.k.noalias() = - llt_.solve(kkt_residual.lu);
  assert(!lqr_policy.K.hasNaN());
  assert(!lqr_policy.k.hasNaN());
  backward_recursion_.factorizeRiccatiFactorization(riccati_next, kkt_matrix, 
                                                    kkt_residual, lqr_policy,
                                                    riccati);
}


inline void RiccatiFactorizer::backwardRiccatiRecursion(
    const SplitRiccatiFactorization& riccati_next, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual, 
    const SplitSwitchingConstraintJacobian& switch_jacobian,
    const SplitSwitchingConstraintResidual& switch_residual, 
    SplitRiccatiFactorization& riccati,
    SplitConstrainedRiccatiFactorization& c_riccati, LQRPolicy& lqr_policy) {
  backward_recursion_.factorizeKKTMatrix(riccati_next, kkt_matrix, kkt_residual);
  // Schur complement
  llt_.compute(kkt_matrix.Quu);
  assert(llt_.info() == Eigen::Success);
  c_riccati.setImpulseStatus(switch_jacobian.dimi());
  c_riccati.Ginv.noalias() = llt_.solve(Eigen::MatrixXd::Identity(dimu_, dimu_));
  c_riccati.DGinv().transpose().noalias() = llt_.solve(switch_jacobian.Phiu().transpose());
  c_riccati.S().noalias() = c_riccati.DGinv() * switch_jacobian.Phiu().transpose();
  llt_s_.compute(c_riccati.S());
  assert(llt_s_.info() == Eigen::Success);
  c_riccati.SinvDGinv().noalias() = llt_s_.solve(c_riccati.DGinv());
  c_riccati.Ginv.noalias() -= c_riccati.SinvDGinv().transpose() * c_riccati.DGinv();
  lqr_policy.K.noalias()  = - c_riccati.Ginv * kkt_matrix.Qxu.transpose();
  lqr_policy.K.noalias() -= c_riccati.SinvDGinv().transpose() * switch_jacobian.Phix();
  lqr_policy.k.noalias()  = - c_riccati.Ginv * kkt_residual.lu;
  lqr_policy.k.noalias() -= c_riccati.SinvDGinv().transpose() * switch_residual.P();
  c_riccati.M().noalias()  = llt_s_.solve(switch_jacobian.Phix());
  c_riccati.M().noalias() -= c_riccati.SinvDGinv() * kkt_matrix.Qxu.transpose();
  c_riccati.m().noalias()  = llt_s_.solve(switch_residual.P());
  c_riccati.m().noalias() -= c_riccati.SinvDGinv() * kkt_residual.lu;
  assert(!lqr_policy.K.hasNaN());
  assert(!lqr_policy.k.hasNaN());
  assert(!c_riccati.M().hasNaN());
  assert(!c_riccati.m().hasNaN());
  backward_recursion_.factorizeRiccatiFactorization(riccati_next, kkt_matrix, 
                                                    kkt_residual, lqr_policy,
                                                    riccati);
  c_riccati.DtM.noalias()   = switch_jacobian.Phiu().transpose() * c_riccati.M();
  c_riccati.KtDtM.noalias() = lqr_policy.K.transpose() * c_riccati.DtM;
  riccati.P.noalias() -= c_riccati.KtDtM;
  riccati.P.noalias() -= c_riccati.KtDtM.transpose();
  riccati.s.noalias() -= switch_jacobian.Phix().transpose() * c_riccati.m();
}


inline void RiccatiFactorizer::backwardRiccatiRecursion(
    const SplitRiccatiFactorization& riccati_next, 
    ImpulseSplitKKTMatrix& kkt_matrix, ImpulseSplitKKTResidual& kkt_residual, 
    SplitRiccatiFactorization& riccati) {
  backward_recursion_.factorizeKKTMatrix(riccati_next, kkt_matrix);
  backward_recursion_.factorizeRiccatiFactorization(riccati_next, kkt_matrix, 
                                                    kkt_residual, riccati);
}


template <typename SplitDirectionType>
inline void RiccatiFactorizer::forwardRiccatiRecursion(
    const SplitKKTMatrix& kkt_matrix, const SplitKKTResidual& kkt_residual, 
    const LQRPolicy& lqr_policy, SplitDirection& d, 
    SplitDirectionType& d_next) const {
  d.du.noalias()  = lqr_policy.K * d.dx;
  d.du.noalias() += lqr_policy.k;
  d_next.dx = kkt_residual.Fx;
  d_next.dx.noalias() += kkt_matrix.Fxx * d.dx;
  d_next.dv().noalias() += kkt_matrix.Fvu * d.du;
}


inline void RiccatiFactorizer::forwardRiccatiRecursion(
    const ImpulseSplitKKTMatrix& kkt_matrix, 
    const ImpulseSplitKKTResidual& kkt_residual, 
    const ImpulseSplitDirection& d, SplitDirection& d_next) const {
  d_next.dx = kkt_residual.Fx;
  d_next.dx.noalias() += kkt_matrix.Fxx * d.dx;
}


template <typename SplitDirectionType>
inline void RiccatiFactorizer::computeCostateDirection(
    const SplitRiccatiFactorization& riccati, SplitDirectionType& d) {
  d.dlmdgmm.noalias() = riccati.P * d.dx - riccati.s;
}


inline void RiccatiFactorizer::computeLagrangeMultiplierDirection(
    const SplitConstrainedRiccatiFactorization& c_riccati, 
    SplitDirection& d) {
  d.dxi().noalias()  = c_riccati.M() * d.dx;
  d.dxi().noalias() += c_riccati.m();
}

} // namespace idocp

#endif // IDOCP_RICCATI_FACTORIZER_HXX_ 