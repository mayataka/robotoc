#include "robotoc/riccati/backward_riccati_recursion_factorizer.hpp"


namespace robotoc {

BackwardRiccatiRecursionFactorizer::BackwardRiccatiRecursionFactorizer(
    const Robot& robot) 
  : dimv_(robot.dimv()),
    dimu_(robot.dimu()),
    AtP_(MatrixXdRowMajor::Zero(2*robot.dimv(), 2*robot.dimv())),
    BtP_(MatrixXdRowMajor::Zero(robot.dimu(), 2*robot.dimv())),
    GK_(Eigen::MatrixXd::Zero(robot.dimu(), 2*robot.dimv())), 
    Pf_(Eigen::VectorXd::Zero(2*robot.dimv())) {
}


BackwardRiccatiRecursionFactorizer::BackwardRiccatiRecursionFactorizer() 
  : dimv_(0),
    dimu_(0),
    AtP_(),
    BtP_(),
    GK_(),
    Pf_() {
}


BackwardRiccatiRecursionFactorizer::~BackwardRiccatiRecursionFactorizer() {
}


void BackwardRiccatiRecursionFactorizer::factorizeKKTMatrix(
    const SplitRiccatiFactorization& riccati_next, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) {
  AtP_.noalias() = kkt_matrix.Fxx.transpose() * riccati_next.P;
  BtP_.noalias() = kkt_matrix.Fvu.transpose() * riccati_next.P.bottomRows(dimv_);
  // Factorize F
  kkt_matrix.Qxx.noalias() += AtP_ * kkt_matrix.Fxx;
  // Factorize H
  kkt_matrix.Qxu.noalias() += AtP_.rightCols(dimv_) * kkt_matrix.Fvu;
  // Factorize G
  kkt_matrix.Quu.noalias() += BtP_.rightCols(dimv_) * kkt_matrix.Fvu;
  // Factorize vector term
  kkt_residual.lu.noalias() += BtP_ * kkt_residual.Fx;
  kkt_residual.lu.noalias() -= kkt_matrix.Fvu.transpose() * riccati_next.sv();
}


void BackwardRiccatiRecursionFactorizer::factorizeHamiltonian(
    const SplitRiccatiFactorization& riccati_next, 
    const SplitKKTMatrix& kkt_matrix, SplitRiccatiFactorization& riccati,
    const bool has_next_sto_phase) const {
  riccati.psi_x.noalias() = AtP_ * kkt_matrix.fx;
  riccati.psi_u.noalias() = BtP_ * kkt_matrix.fx;
  riccati.psi_x.noalias() += kkt_matrix.hx;
  riccati.psi_u.noalias() += kkt_matrix.hu;
  riccati.psi_x.noalias() += kkt_matrix.Fxx.transpose() * riccati_next.Psi;
  riccati.psi_u.noalias() += kkt_matrix.Fvu.transpose() * riccati_next.Psi.tail(dimv_);
  if (has_next_sto_phase) {
    riccati.phi_x.noalias() = kkt_matrix.Fxx.transpose() * riccati_next.Phi;
    riccati.phi_u.noalias() = kkt_matrix.Fvu.transpose() * riccati_next.Phi.tail(dimv_);
  }
  else {
    riccati.phi_x.setZero();
    riccati.phi_u.setZero();
  }
}


void BackwardRiccatiRecursionFactorizer::factorizeKKTMatrix(
    const SplitRiccatiFactorization& riccati_next, 
    SplitKKTMatrix& kkt_matrix) {
  AtP_.noalias() = kkt_matrix.Fxx.transpose() * riccati_next.P;
  // Factorize F
  kkt_matrix.Qxx.noalias() += AtP_ * kkt_matrix.Fxx;
}


void BackwardRiccatiRecursionFactorizer::factorizeRiccatiFactorization(
    const SplitRiccatiFactorization& riccati_next, SplitKKTMatrix& kkt_matrix, 
    const SplitKKTResidual& kkt_residual, const LQRPolicy& lqr_policy, 
    SplitRiccatiFactorization& riccati) {
  GK_.noalias() = kkt_matrix.Quu * lqr_policy.K; 
  kkt_matrix.Qxx.noalias() -= lqr_policy.K.transpose() * GK_;
  // Riccati factorization matrix with preserving the symmetry
  riccati.P = 0.5 * (kkt_matrix.Qxx + kkt_matrix.Qxx.transpose());
  // Riccati factorization vector
  riccati.s.noalias()  = kkt_matrix.Fxx.transpose() * riccati_next.s;
  riccati.s.noalias() -= AtP_ * kkt_residual.Fx;
  riccati.s.noalias() -= kkt_residual.lx;
  riccati.s.noalias() -= kkt_matrix.Qxu * lqr_policy.k;
}


void BackwardRiccatiRecursionFactorizer::factorizeSTOFactorization(
    const SplitRiccatiFactorization& riccati_next, 
    const SplitKKTMatrix& kkt_matrix, const SplitKKTResidual& kkt_residual, 
    const LQRPolicy& lqr_policy, SplitRiccatiFactorization& riccati,
    const bool has_next_sto_phase) {
  // Qtx
  riccati.Psi            = riccati.psi_x;
  riccati.Psi.noalias() += lqr_policy.K.transpose() * riccati.psi_u;
  if (has_next_sto_phase) {
    riccati.Phi.noalias()  = riccati.phi_x;
    riccati.Phi.noalias() += lqr_policy.K.transpose() * riccati.phi_u;
  }
  else {
    riccati.Phi.setZero();
  }
  // Qtt
  Pf_.noalias() = riccati_next.P * kkt_matrix.fx;
  riccati.xi    = kkt_matrix.fx.dot(Pf_);
  riccati.xi   += kkt_matrix.Qtt;
  riccati.xi   += 2 * riccati_next.Psi.dot(kkt_matrix.fx);
  riccati.xi   += lqr_policy.T.dot(riccati.psi_u);
  riccati.xi   += riccati_next.xi;
  if (has_next_sto_phase) {
    riccati.chi  = kkt_matrix.Qtt_prev;
    riccati.chi += riccati_next.Phi.dot(kkt_matrix.fx);
    riccati.chi += lqr_policy.T.dot(riccati.phi_u);
    riccati.chi += riccati_next.chi;
    riccati.rho  = lqr_policy.W.dot(riccati.phi_u);
    riccati.rho += riccati_next.rho;
  }
  else {
    riccati.chi = 0.0;
    riccati.rho = 0.0;
  }
  // h
  Pf_.noalias() = riccati_next.P * kkt_residual.Fx - riccati_next.s;
  riccati.eta   = kkt_matrix.fx.dot(Pf_);
  riccati.eta  += kkt_residual.h;
  riccati.eta  += riccati_next.Psi.dot(kkt_residual.Fx);
  riccati.eta  += riccati.psi_u.dot(lqr_policy.k);
  riccati.eta  += riccati_next.eta;
  if (has_next_sto_phase) {
    riccati.iota  = riccati_next.Phi.dot(kkt_residual.Fx);
    riccati.iota += riccati.phi_u.dot(lqr_policy.k);
    riccati.iota += riccati_next.iota;
  }
  else {
    riccati.iota = 0.0;
  }
}


void BackwardRiccatiRecursionFactorizer::factorizeRiccatiFactorization(
    const SplitRiccatiFactorization& riccati_next, 
    const SplitKKTMatrix& kkt_matrix, 
    const SplitKKTResidual& kkt_residual, 
    SplitRiccatiFactorization& riccati) {
  // Riccati factorization matrix with preserving the symmetry
  riccati.P = 0.5 * (kkt_matrix.Qxx + kkt_matrix.Qxx.transpose());
  // Riccati factorization vector
  riccati.s.noalias()  = kkt_matrix.Fxx.transpose() * riccati_next.s;
  riccati.s.noalias() -= AtP_ * kkt_residual.Fx;
  riccati.s.noalias() -= kkt_residual.lx;
}


void BackwardRiccatiRecursionFactorizer::factorizeSTOFactorization(
    const SplitRiccatiFactorization& riccati_next, 
    const SplitKKTMatrix& kkt_matrix, 
    const SplitKKTResidual& kkt_residual, 
    SplitRiccatiFactorization& riccati) {
  // Qtx
  riccati.Psi.setZero();
  riccati.Phi.noalias() = kkt_matrix.Fxx.transpose() * riccati_next.Phi;
  riccati.xi    = 0.0; 
  riccati.chi   = 0.0;
  riccati.rho   = riccati_next.rho;
  riccati.eta   = 0.0;
  riccati.iota  = riccati_next.iota;
  riccati.iota += riccati_next.Phi.dot(kkt_residual.Fx);
}

} // namespace robotoc