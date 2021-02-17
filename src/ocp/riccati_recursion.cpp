#include "idocp/ocp/riccati_recursion.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>

namespace idocp {

RiccatiRecursion::RiccatiRecursion(const Robot& robot, const int N, 
                                   const int max_num_impulse)
  : N_(N),
    factorizer_(robot, N, max_num_impulse) {
  try {
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
    if (max_num_impulse < 0) {
      throw std::out_of_range("invalid value: max_num_impulse must be non-negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


RiccatiRecursion::RiccatiRecursion()
  : N_(0),
    factorizer_() {
}


RiccatiRecursion::~RiccatiRecursion() {
}


void RiccatiRecursion::backwardRiccatiRecursion(
    const OCPDiscretizer& ocp_discretizer, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual, const StateConstraintJacobian& jac,
    RiccatiFactorization& factorization) {
  factorization[N_].Pqq = kkt_matrix[N_].Qqq();
  factorization[N_].Pvv = kkt_matrix[N_].Qvv();
  factorization[N_].sq = - kkt_residual[N_].lq();
  factorization[N_].sv = - kkt_residual[N_].lv();
  for (int i=N_-1; i>=0; --i) {
    if (ocp_discretizer.isTimeStageBeforeImpulse(i)) {
      const int impulse_index = ocp_discretizer.impulseIndexAfterTimeStage(i);
      factorizer_.aux[impulse_index].backwardRiccatiRecursion(
          factorization[i+1], ocp_discretizer.dtau_aux(impulse_index), 
          kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index], 
          factorization.aux[impulse_index]);
      factorizer_.impulse[impulse_index].backwardRiccatiRecursion(
          factorization.aux[impulse_index],  kkt_matrix.impulse[impulse_index], 
          kkt_residual.impulse[impulse_index], 
          factorization.impulse[impulse_index]);
      factorizer_[i].backwardRiccatiRecursion(
          factorization.impulse[impulse_index], ocp_discretizer.dtau(i), 
          kkt_matrix[i], kkt_residual[i], factorization[i]);
    }
    else if (ocp_discretizer.isTimeStageBeforeLift(i)) {
      const int lift_index = ocp_discretizer.liftIndexAfterTimeStage(i);
      if (ocp_discretizer.isTimeStageBeforeImpulse(i+1)) {
        const int impulse_index = ocp_discretizer.impulseIndexAfterTimeStage(i+1);
        factorizer_.lift[lift_index].backwardRiccatiRecursion(
            factorization[i+1], ocp_discretizer.dtau_lift(lift_index), 
            kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index], 
            jac[impulse_index], factorization.lift[lift_index],
            factorization.constraint[impulse_index]);
      }
      else {
        factorizer_.lift[lift_index].backwardRiccatiRecursion(
            factorization[i+1], ocp_discretizer.dtau_lift(lift_index), 
            kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index], 
            factorization.lift[lift_index]);
      }
      factorizer_[i].backwardRiccatiRecursion(
          factorization.lift[lift_index], ocp_discretizer.dtau(i), 
          kkt_matrix[i], kkt_residual[i], factorization[i]);
    }
    else {
      if (ocp_discretizer.isTimeStageBeforeImpulse(i+1)) {
        const int impulse_index = ocp_discretizer.impulseIndexAfterTimeStage(i+1);
        factorizer_[i].backwardRiccatiRecursion(
            factorization[i+1], ocp_discretizer.dtau(i), kkt_matrix[i], 
            kkt_residual[i], jac[impulse_index], factorization[i], 
            factorization.constraint[impulse_index]);
      }
      else {
        factorizer_[i].backwardRiccatiRecursion(
            factorization[i+1], ocp_discretizer.dtau(i), 
            kkt_matrix[i], kkt_residual[i], factorization[i]);
      }
    }
  }
}


void RiccatiRecursion::forwardRiccatiRecursion(
    const OCPDiscretizer& ocp_discretizer, const KKTMatrix& kkt_matrix, 
    const KKTResidual& kkt_residual, Direction& d) const {
  for (int i=0; i<N_; ++i) {
    if (ocp_discretizer.isTimeStageBeforeImpulse(i)) {
      const int impulse_index = ocp_discretizer.impulseIndexAfterTimeStage(i);
      factorizer_[i].forwardRiccatiRecursion(kkt_matrix[i], kkt_residual[i],  
                                             ocp_discretizer.dtau(i), d[i],
                                             d.impulse[impulse_index]);
      factorizer_.impulse[impulse_index].forwardRiccatiRecursion(
          kkt_matrix.impulse[impulse_index], kkt_residual.impulse[impulse_index],
          d.impulse[impulse_index], d.aux[impulse_index]);
      factorizer_.aux[impulse_index].forwardRiccatiRecursion(
          kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index],
          ocp_discretizer.dtau_aux(impulse_index), d.aux[impulse_index], d[i+1]);
    }
    else if (ocp_discretizer.isTimeStageBeforeLift(i)) {
      const int lift_index = ocp_discretizer.liftIndexAfterTimeStage(i);
      factorizer_[i].forwardRiccatiRecursion(kkt_matrix[i], kkt_residual[i],  
                                             ocp_discretizer.dtau(i), d[i],
                                             d.lift[lift_index]);
      factorizer_.lift[lift_index].forwardRiccatiRecursion(
          kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index],
          ocp_discretizer.dtau_lift(lift_index), d.lift[lift_index], d[i+1]);
    }
    else {
      factorizer_[i].forwardRiccatiRecursion(kkt_matrix[i], kkt_residual[i],  
                                             ocp_discretizer.dtau(i), 
                                             d[i], d[i+1]);
    }
  }
}

} // namespace idocp