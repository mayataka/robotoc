#include "robotoc/riccati/unconstr_riccati_recursion.hpp"

#include <omp.h>
#include <stdexcept>
#include <iostream>
#include <cassert>

namespace robotoc {

UnconstrRiccatiRecursion::UnconstrRiccatiRecursion(const UnconstrOCP& ocp)
  : N_(ocp.N()),
    dt_(ocp.T()/ocp.N()),
    factorizer_(ocp.robot()),
    lqr_policy_(ocp.N(), LQRPolicy(ocp.robot())) {
}


UnconstrRiccatiRecursion::UnconstrRiccatiRecursion()
  : N_(0),
    dt_(0),
    factorizer_(),
    lqr_policy_() {
}


UnconstrRiccatiRecursion::~UnconstrRiccatiRecursion() {
}


void UnconstrRiccatiRecursion::backwardRiccatiRecursion(
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual, 
    UnconstrRiccatiFactorization& factorization) {
  factorization[N_].P = kkt_matrix[N_].Qxx;
  factorization[N_].s = - kkt_residual[N_].lx;
  for (int i=N_-1; i>=0; --i) {
    factorizer_.backwardRiccatiRecursion(factorization[i+1], dt_, kkt_matrix[i], 
                                         kkt_residual[i], factorization[i], 
                                         lqr_policy_[i]);
  }
}


void UnconstrRiccatiRecursion::forwardRiccatiRecursion( 
    const KKTResidual& kkt_residual, Direction& d) const {
  for (int i=0; i<N_; ++i) {
    factorizer_.forwardRiccatiRecursion(kkt_residual[i], dt_, lqr_policy_[i], 
                                        d[i], d[i+1]);
  }
}


void UnconstrRiccatiRecursion::getStateFeedbackGain(
    const int time_stage, Eigen::MatrixXd& da_dq, 
    Eigen::MatrixXd& da_dv) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  da_dq = lqr_policy_[time_stage].Kq();
  da_dv = lqr_policy_[time_stage].Kv();
}

} // namespace robotoc