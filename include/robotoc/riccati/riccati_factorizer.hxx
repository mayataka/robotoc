#ifndef ROBOTOC_RICCATI_FACTORIZER_HXX_ 
#define ROBOTOC_RICCATI_FACTORIZER_HXX_

#include "robotoc/riccati/riccati_factorizer.hpp"

#include <cassert>

namespace robotoc {

inline RiccatiFactorizer::RiccatiFactorizer(const Robot& robot, 
                                            const double max_dts0) 
  : has_floating_base_(robot.hasFloatingBase()),
    dimv_(robot.dimv()),
    dimu_(robot.dimu()),
    max_dts0_(max_dts0),
    llt_(robot.dimu()),
    llt_s_(),
    backward_recursion_(robot) {
}


inline RiccatiFactorizer::RiccatiFactorizer() 
  : has_floating_base_(false),
    dimv_(0),
    dimu_(0),
    max_dts0_(1.0),
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
    const SplitRiccatiFactorization& riccati_next,  SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual, SplitRiccatiFactorization& riccati, 
    LQRPolicy& lqr_policy, const bool sto, const bool has_next_sto_phase) {
  backwardRiccatiRecursion(riccati_next, kkt_matrix, kkt_residual,
                           riccati, lqr_policy);
  if (sto) {
    backward_recursion_.factorizeHamiltonian(riccati_next, kkt_matrix, riccati,
                                             has_next_sto_phase);
    lqr_policy.T.noalias() = - llt_.solve(riccati.psi_u);
    if (has_next_sto_phase) {
      lqr_policy.W.noalias() = - llt_.solve(riccati.phi_u);
    }
    else {
      lqr_policy.W.setZero();
    }
    backward_recursion_.factorizeSTOFactorization(riccati_next, kkt_matrix, 
                                                  kkt_residual, lqr_policy, 
                                                  riccati, has_next_sto_phase);
  }
  else {
    riccati.Psi.setZero();
    riccati.xi = 0.;
    riccati.chi = 0.;
    riccati.eta = 0.;
  }
}


inline void RiccatiFactorizer::backwardRiccatiRecursionPhaseTransition(
    const SplitRiccatiFactorization& riccati, 
    SplitRiccatiFactorization& riccati_m, STOPolicy& sto_policy, 
    const bool has_next_sto_phase) const {
  riccati_m.P = riccati.P;
  riccati_m.s = riccati.s;
  riccati_m.Psi.setZero();
  riccati_m.Phi = riccati.Psi;
  riccati_m.xi = 0.0;
  riccati_m.chi = 0.0;
  riccati_m.rho = riccati.xi;
  riccati_m.eta = 0.0;
  riccati_m.iota = riccati.eta;
  if (has_next_sto_phase) {
    double sgm = riccati.xi - 2.0 * riccati.chi + riccati.rho;
    if ((sgm*max_dts0_) < std::abs(riccati.eta-riccati.iota) || sgm < keps_) {
      std::cout << "sgm reg ! sgm = " << sgm << std::endl;
      std::cout << "sgm * max_dts0_ = " << sgm*max_dts0_ << std::endl;
      std::cout << "std::abs(riccati.eta-riccati.iota) = " << std::abs(riccati.eta-riccati.iota) << std::endl;
      sgm = std::abs(sgm) + std::abs(riccati.eta-riccati.iota) / max_dts0_;
    }
    sto_policy.dtsdx  = - (1.0/sgm) * (riccati.Psi-riccati.Phi);
    sto_policy.dtsdts =   (1.0/sgm) * (riccati.xi-riccati.chi);
    sto_policy.dts0   = - (1.0/sgm) * (riccati.eta-riccati.iota);
    riccati_m.s.noalias()   
        += (1.0/sgm) * (riccati.Psi-riccati.Phi) * (riccati.eta-riccati.iota);
    riccati_m.Phi.noalias() 
        -= (1.0/sgm) * (riccati.Psi-riccati.Phi) * (riccati.xi-riccati.chi);
    riccati_m.rho
        = riccati.xi - (1.0/sgm) * (riccati.xi-riccati.chi) * (riccati.xi-riccati.chi);
    riccati_m.iota
        = riccati.eta - (1.0/sgm) * (riccati.xi-riccati.chi)  * (riccati.eta-riccati.iota);
  }
}


inline void RiccatiFactorizer::backwardRiccatiRecursion(
    const SplitRiccatiFactorization& riccati_next, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual, 
    const SwitchingConstraintJacobian& sc_jacobian,
    const SwitchingConstraintResidual& sc_residual, 
    SplitRiccatiFactorization& riccati,
    SplitConstrainedRiccatiFactorization& c_riccati, LQRPolicy& lqr_policy) {
  backward_recursion_.factorizeKKTMatrix(riccati_next, kkt_matrix, kkt_residual);
  // Schur complement
  llt_.compute(kkt_matrix.Quu);
  assert(llt_.info() == Eigen::Success);
  c_riccati.setImpulseStatus(sc_jacobian.dimi());
  c_riccati.Ginv.noalias() = llt_.solve(Eigen::MatrixXd::Identity(dimu_, dimu_));
  c_riccati.DGinv().transpose().noalias() = llt_.solve(sc_jacobian.Phiu().transpose());
  c_riccati.S().noalias() = c_riccati.DGinv() * sc_jacobian.Phiu().transpose();
  llt_s_.compute(c_riccati.S());
  assert(llt_s_.info() == Eigen::Success);
  c_riccati.SinvDGinv().noalias() = llt_s_.solve(c_riccati.DGinv());
  c_riccati.Ginv.noalias() -= c_riccati.SinvDGinv().transpose() * c_riccati.DGinv();
  lqr_policy.K.noalias()  = - c_riccati.Ginv * kkt_matrix.Qxu.transpose();
  lqr_policy.K.noalias() -= c_riccati.SinvDGinv().transpose() * sc_jacobian.Phix();
  lqr_policy.k.noalias()  = - c_riccati.Ginv * kkt_residual.lu;
  lqr_policy.k.noalias() -= c_riccati.SinvDGinv().transpose() * sc_residual.P();
  c_riccati.M().noalias()  = llt_s_.solve(sc_jacobian.Phix());
  c_riccati.M().noalias() -= c_riccati.SinvDGinv() * kkt_matrix.Qxu.transpose();
  c_riccati.m().noalias()  = llt_s_.solve(sc_residual.P());
  c_riccati.m().noalias() -= c_riccati.SinvDGinv() * kkt_residual.lu;
  assert(!lqr_policy.K.hasNaN());
  assert(!lqr_policy.k.hasNaN());
  assert(!c_riccati.M().hasNaN());
  assert(!c_riccati.m().hasNaN());
  backward_recursion_.factorizeRiccatiFactorization(riccati_next, kkt_matrix, 
                                                    kkt_residual, lqr_policy,
                                                    riccati);
  c_riccati.DtM.noalias()   = sc_jacobian.Phiu().transpose() * c_riccati.M();
  c_riccati.KtDtM.noalias() = lqr_policy.K.transpose() * c_riccati.DtM;
  riccati.P.noalias() -= c_riccati.KtDtM;
  riccati.P.noalias() -= c_riccati.KtDtM.transpose();
  riccati.s.noalias() -= sc_jacobian.Phix().transpose() * c_riccati.m();
}


inline void RiccatiFactorizer::backwardRiccatiRecursion(
    const SplitRiccatiFactorization& riccati_next, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual, 
    const SwitchingConstraintJacobian& sc_jacobian,
    const SwitchingConstraintResidual& sc_residual, 
    SplitRiccatiFactorization& riccati,
    SplitConstrainedRiccatiFactorization& c_riccati, LQRPolicy& lqr_policy,
    const bool sto, const bool has_next_sto_phase) {
  backwardRiccatiRecursion(riccati_next, kkt_matrix, kkt_residual, sc_jacobian, 
                           sc_residual, riccati, c_riccati, lqr_policy);
  if (sto) {
    backward_recursion_.factorizeHamiltonian(riccati_next, kkt_matrix, riccati,
                                             has_next_sto_phase);
    lqr_policy.T.noalias()  = - c_riccati.Ginv * riccati.psi_u;
    lqr_policy.T.noalias() -= c_riccati.SinvDGinv().transpose() * sc_jacobian.Phit();
    if (has_next_sto_phase) {
      lqr_policy.W.noalias()  = - c_riccati.Ginv * riccati.phi_u;
    }
    else {
      lqr_policy.W.setZero();
    }
    c_riccati.mt().noalias()  = llt_s_.solve(sc_jacobian.Phit());
    c_riccati.mt().noalias() -= c_riccati.SinvDGinv() * riccati.psi_u;
    if (has_next_sto_phase) {
      c_riccati.mt_next().noalias() = - c_riccati.SinvDGinv() * riccati.phi_u;
    }
    else {
      c_riccati.mt_next().setZero();
    }
    backward_recursion_.factorizeSTOFactorization(riccati_next, kkt_matrix, 
                                                  kkt_residual, lqr_policy, 
                                                  riccati, has_next_sto_phase);
    riccati.Psi.noalias() += c_riccati.M().transpose() * sc_jacobian.Phit();
    riccati.xi += c_riccati.mt().dot(sc_jacobian.Phit());
    if (has_next_sto_phase) {
      riccati.chi += c_riccati.mt_next().dot(sc_jacobian.Phit());
    }
    else {
      riccati.chi = 0.0;
    }
    riccati.eta += c_riccati.m().dot(sc_jacobian.Phit());
  }
}


inline void RiccatiFactorizer::backwardRiccatiRecursion(
    const SplitRiccatiFactorization& riccati_next, 
    ImpulseSplitKKTMatrix& kkt_matrix, ImpulseSplitKKTResidual& kkt_residual, 
    SplitRiccatiFactorization& riccati) {
  backward_recursion_.factorizeKKTMatrix(riccati_next, kkt_matrix);
  backward_recursion_.factorizeRiccatiFactorization(riccati_next, kkt_matrix, 
                                                    kkt_residual, riccati);
}


inline void RiccatiFactorizer::backwardRiccatiRecursion(
    const SplitRiccatiFactorization& riccati_next, 
    ImpulseSplitKKTMatrix& kkt_matrix, ImpulseSplitKKTResidual& kkt_residual, 
    SplitRiccatiFactorization& riccati, const bool sto) {
  backwardRiccatiRecursion(riccati_next, kkt_matrix, kkt_residual, riccati);
  if (sto) {
    backward_recursion_.factorizeSTOFactorization(riccati_next, kkt_matrix, 
                                                  kkt_residual, riccati);
  }
}


inline void RiccatiFactorizer::computeSwitchingTimeDirection(
    const STOPolicy& sto_policy, SplitDirection& d, 
    const bool has_prev_sto_phase) {
  d.dts_next = sto_policy.dtsdx.dot(d.dx) + sto_policy.dts0;
  if (has_prev_sto_phase) {
    d.dts_next += sto_policy.dtsdts * d.dts;
  }
}


template <typename SplitDirectionType>
inline void RiccatiFactorizer::forwardRiccatiRecursion(
    const SplitKKTMatrix& kkt_matrix, const SplitKKTResidual& kkt_residual, 
    const LQRPolicy& lqr_policy, SplitDirection& d, SplitDirectionType& d_next, 
    const bool sto, const bool has_next_sto_phase) const {
  d.du.noalias()  = lqr_policy.K * d.dx;
  d.du.noalias() += lqr_policy.k;
  if (sto) {
    d.du.noalias() += lqr_policy.T * (d.dts_next-d.dts);
    if (has_next_sto_phase) {
      d.du.noalias() -= lqr_policy.W * d.dts_next;
    }
  }
  d_next.dx = kkt_residual.Fx;
  d_next.dx.noalias()   += kkt_matrix.Fxx * d.dx;
  d_next.dv().noalias() += kkt_matrix.Fvu * d.du;
  if (sto) {
    d_next.dx.noalias() += kkt_matrix.fx * (d.dts_next-d.dts);
  }
  d_next.dts = d.dts;
  d_next.dts_next = d.dts_next;
}


inline void RiccatiFactorizer::forwardRiccatiRecursion(
    const ImpulseSplitKKTMatrix& kkt_matrix, 
    const ImpulseSplitKKTResidual& kkt_residual, 
    const ImpulseSplitDirection& d, SplitDirection& d_next) const {
  d_next.dx = kkt_residual.Fx;
  d_next.dx.noalias() += kkt_matrix.Fxx * d.dx;
  d_next.dts = d.dts;
  d_next.dts_next = d.dts_next;
}


inline void RiccatiFactorizer::computeCostateDirection(
    const SplitRiccatiFactorization& riccati, SplitDirection& d,
    const bool sto, const bool has_next_sto_phase) {
  d.dlmdgmm.noalias() = riccati.P * d.dx - riccati.s;
  if (sto) {
    d.dlmdgmm.noalias() += riccati.Psi * (d.dts_next-d.dts);
    if (has_next_sto_phase) {
      d.dlmdgmm.noalias() -= riccati.Phi * d.dts_next;
    }
  }
}


inline void RiccatiFactorizer::computeCostateDirection(
    const SplitRiccatiFactorization& riccati, ImpulseSplitDirection& d,
    const bool sto) {
  d.dlmdgmm.noalias() = riccati.P * d.dx - riccati.s;
  if (sto) {
    d.dlmdgmm.noalias() -= riccati.Phi * d.dts_next;
  }
}


inline void RiccatiFactorizer::computeLagrangeMultiplierDirection(
    const SplitConstrainedRiccatiFactorization& c_riccati, 
    SplitDirection& d, const bool sto, const bool has_next_sto_phase) {
  d.dxi().noalias()  = c_riccati.M() * d.dx;
  d.dxi().noalias() += c_riccati.m();
  if (sto) {
    d.dxi().noalias() += c_riccati.mt() * (d.dts_next-d.dts);
    if (has_next_sto_phase) {
      d.dxi().noalias() -= c_riccati.mt_next() * d.dts_next;
    }
  }
}

} // namespace robotoc

#endif // ROBOTOC_RICCATI_FACTORIZER_HXX_ 