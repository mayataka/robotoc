#include "idocp/riccati/unconstr_riccati_recursion.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>

namespace idocp {

UnconstrRiccatiRecursion::UnconstrRiccatiRecursion(const Robot& robot, 
                                                   const double T, const int N)
  : N_(N),
    T_(T),
    dt_(T/N),
    factorizer_(robot),
    lqr_policy_(N, LQRPolicy(robot)) {
  try {
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


UnconstrRiccatiRecursion::UnconstrRiccatiRecursion()
  : N_(0),
    T_(0),
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


void UnconstrRiccatiRecursion::getStateFeedbackGain(const int time_stage, 
                                                    Eigen::MatrixXd& Kq, 
                                                    Eigen::MatrixXd& Kv) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  Kq = lqr_policy_[time_stage].Kq();
  Kv = lqr_policy_[time_stage].Kv();
}

} // namespace idocp