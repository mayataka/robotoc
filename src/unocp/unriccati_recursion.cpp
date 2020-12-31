#include "idocp/unocp/unriccati_recursion.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>

namespace idocp {

UnRiccatiRecursion::UnRiccatiRecursion(const Robot& robot, const double T, 
                                       const int N)
  : N_(N),
    T_(T),
    dtau_(T/N),
    factorizer_(N, SplitUnRiccatiFactorizer(robot)) {
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


UnRiccatiRecursion::UnRiccatiRecursion()
  : N_(0),
    T_(0),
    dtau_(0),
    factorizer_() {
}


UnRiccatiRecursion::~UnRiccatiRecursion() {
}


void UnRiccatiRecursion::backwardRiccatiRecursionTerminal(
    const SplitKKTMatrix& terminal_kkt_matrix, 
    const SplitKKTResidual& terminal_kkt_residual,
    UnRiccatiFactorization& riccati_factorization) const {
  riccati_factorization[N_].Pqq = terminal_kkt_matrix.Qqq();
  riccati_factorization[N_].Pvv = terminal_kkt_matrix.Qvv();
  riccati_factorization[N_].sq = - terminal_kkt_residual.lq();
  riccati_factorization[N_].sv = - terminal_kkt_residual.lv();
}


void UnRiccatiRecursion::backwardRiccatiRecursion(
    UnKKTMatrix& unkkt_matrix, UnKKTResidual& unkkt_residual, 
    UnRiccatiFactorization& riccati_factorization) {
  for (int i=N_-1; i>=0; --i) {
    factorizer_[i].backwardRiccatiRecursion(
        riccati_factorization[i+1], dtau_, 
        unkkt_matrix[i], unkkt_residual[i], riccati_factorization[i]);
  }
}

void UnRiccatiRecursion::forwardRiccatiRecursion( 
    const UnKKTResidual& unkkt_residual, UnDirection& d) const {
  for (int i=0; i<N_; ++i) {
    factorizer_[i].forwardRiccatiRecursion(unkkt_residual[i], d[i], dtau_, 
                                           d[i+1]);
  }
}

} // namespace idocp