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
    RiccatiFactorization& riccati_factorization) {
  riccati_factorization[N_].Pqq = kkt_matrix[N_].Qqq();
  riccati_factorization[N_].Pvv = kkt_matrix[N_].Qvv();
  riccati_factorization[N_].sq = - kkt_residual[N_].lq();
  riccati_factorization[N_].sv = - kkt_residual[N_].lv();
  for (int i=N_-1; i>=0; --i) {
    factorizer_.backwardRiccatiRecursion(riccati_factorization[i+1], dt_, 
                                         kkt_matrix[i], kkt_residual[i], 
                                         riccati_factorization[i], 
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

} // namespace idocp