#include "robotoc/riccati/riccati_recursion.hpp"

#include <omp.h>
#include <stdexcept>
#include <iostream>
#include <cassert>

namespace robotoc {

RiccatiRecursion::RiccatiRecursion(const Robot& robot, const int N, 
                                   const int max_num_impulse, 
                                   const int nthreads)
  : nthreads_(nthreads),
    N_(N),
    N_all_(N+1),
    factorizer_(robot),
    lqr_policy_(robot, N, max_num_impulse),
    max_primal_step_sizes_(Eigen::VectorXd::Zero(N+1+3*max_num_impulse)), 
    max_dual_step_sizes_(Eigen::VectorXd::Zero(N+1+3*max_num_impulse)) {
  try {
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
    if (max_num_impulse < 0) {
      throw std::out_of_range("invalid value: max_num_impulse must be non-negative!");
    }
    if (nthreads <= 0) {
      throw std::out_of_range("invalid value: nthreads must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


RiccatiRecursion::RiccatiRecursion()
  : nthreads_(0),
    N_(0),
    N_all_(0),
    factorizer_(),
    max_primal_step_sizes_(), 
    max_dual_step_sizes_() {
}


RiccatiRecursion::~RiccatiRecursion() {
}


void RiccatiRecursion::backwardRiccatiRecursion(
    const OCP& ocp, KKTMatrix& kkt_matrix, KKTResidual& kkt_residual, 
    RiccatiFactorization& factorization) {
  const int N = ocp.discrete().N();
  factorization[N].P = kkt_matrix[N].Qxx;
  factorization[N].s = - kkt_residual[N].lx;
  for (int i=N-1; i>=0; --i) {
    if (ocp.discrete().isTimeStageBeforeImpulse(i)) {
      assert(!ocp.discrete().isTimeStageBeforeImpulse(i+1));
      const int impulse_index = ocp.discrete().impulseIndexAfterTimeStage(i);
      const bool sto = ocp.discrete().isSTOEnabledImpulse(impulse_index);
      factorizer_.backwardRiccatiRecursion(factorization[i+1], 
                                           kkt_matrix.aux[impulse_index], 
                                           kkt_residual.aux[impulse_index], 
                                           factorization.aux[impulse_index], 
                                           lqr_policy_.aux[impulse_index]);
      factorizer_.backwardRiccatiRecursion(factorization.aux[impulse_index], 
                                           kkt_matrix.impulse[impulse_index], 
                                           kkt_residual.impulse[impulse_index], 
                                           factorization.impulse[impulse_index]);
      factorizer_.backwardRiccatiRecursion(factorization.impulse[impulse_index], 
                                           kkt_matrix[i], kkt_residual[i], 
                                           factorization[i], lqr_policy_[i]);
      if (i-1 >= 0) {
        factorizer_.backwardRiccatiRecursion(factorization[i], kkt_matrix[i-1], 
                                             kkt_residual[i-1],
                                             kkt_matrix.switching[impulse_index], 
                                             kkt_residual.switching[impulse_index], 
                                             factorization[i-1], 
                                             factorization.switching[impulse_index], 
                                             lqr_policy_[i-1]);
      }

    }
    else if (ocp.discrete().isTimeStageBeforeLift(i)) {
      assert(!ocp.discrete().isTimeStageBeforeImpulse(i+1));
      const int lift_index = ocp.discrete().liftIndexAfterTimeStage(i);
      const bool sto = ocp.discrete().isSTOEnabledLift(lift_index);
      factorizer_.backwardRiccatiRecursion(factorization[i+1], 
                                            kkt_matrix.lift[lift_index], 
                                            kkt_residual.lift[lift_index], 
                                            factorization.lift[lift_index], 
                                            lqr_policy_.lift[lift_index]);
      factorizer_.backwardRiccatiRecursion(factorization.lift[lift_index], 
                                           kkt_matrix[i], kkt_residual[i], 
                                           factorization[i], lqr_policy_[i]);
    }
    else if (!ocp.discrete().isTimeStageBeforeImpulse(i+1)) {
      factorizer_.backwardRiccatiRecursion(factorization[i+1], 
                                            kkt_matrix[i], kkt_residual[i], 
                                            factorization[i], lqr_policy_[i]);
    }
  }
}


void RiccatiRecursion::forwardRiccatiRecursion(
    const OCP& ocp, const KKTMatrix& kkt_matrix, 
    const KKTResidual& kkt_residual, Direction& d) const {
  const int N = ocp.discrete().N();
  for (int i=0; i<N; ++i) {
    if (ocp.discrete().isTimeStageBeforeImpulse(i)) {
      assert(!ocp.discrete().isTimeStageBeforeImpulse(i+1));
      const int impulse_index = ocp.discrete().impulseIndexAfterTimeStage(i);
      const bool sto = ocp.discrete().isSTOEnabledImpulse(impulse_index);
      if (i-1 >= 0) {
        factorizer_.forwardRiccatiRecursion(kkt_matrix[i-1], kkt_residual[i-1],
                                            lqr_policy_[i-1], d[i-1], d[i]);
      }
      factorizer_.forwardRiccatiRecursion(kkt_matrix[i], kkt_residual[i],  
                                          lqr_policy_[i], d[i], 
                                          d.impulse[impulse_index]);
      factorizer_.forwardRiccatiRecursion(kkt_matrix.impulse[impulse_index], 
                                          kkt_residual.impulse[impulse_index],
                                          d.impulse[impulse_index], d.aux[impulse_index]);
      factorizer_.forwardRiccatiRecursion(kkt_matrix.aux[impulse_index], 
                                          kkt_residual.aux[impulse_index],
                                          lqr_policy_.aux[impulse_index], 
                                          d.aux[impulse_index], d[i+1]);
    }
    else if (ocp.discrete().isTimeStageBeforeLift(i)) {
      assert(!ocp.discrete().isTimeStageBeforeImpulse(i+1));
      const int lift_index = ocp.discrete().liftIndexAfterTimeStage(i);
      const bool sto = ocp.discrete().isSTOEnabledLift(lift_index);
      factorizer_.forwardRiccatiRecursion(kkt_matrix[i], kkt_residual[i],  
                                          lqr_policy_[i], d[i], 
                                          d.lift[lift_index]);
      factorizer_.forwardRiccatiRecursion(kkt_matrix.lift[lift_index], 
                                          kkt_residual.lift[lift_index], 
                                          lqr_policy_.lift[lift_index], 
                                          d.lift[lift_index], d[i+1]);
    }
    else if (!ocp.discrete().isTimeStageBeforeImpulse(i+1)) {
      factorizer_.forwardRiccatiRecursion(kkt_matrix[i], kkt_residual[i],  
                                          lqr_policy_[i], d[i], d[i+1]);
    }
  }
}


void RiccatiRecursion::computeDirection(
    OCP& ocp, const RiccatiFactorization& factorization, 
    const Solution& s, Direction& d) {
  const int N = ocp.discrete().N();
  const int N_impulse = ocp.discrete().N_impulse();
  const int N_lift = ocp.discrete().N_lift();
  const int N_all = N + 1 + 2*N_impulse + N_lift;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N) {
      if (ocp.discrete().isTimeStageBeforeImpulse(i)) {
        const int impulse_index = ocp.discrete().impulseIndexAfterTimeStage(i);
        const bool sto = ocp.discrete().isSTOEnabledImpulse(impulse_index);
        RiccatiFactorizer::computeCostateDirection(factorization[i], d[i], sto);
        ocp[i].expandPrimal(ocp.discrete().dt(i), s[i], d[i], sto);
      }
      else if (ocp.discrete().isTimeStageBeforeLift(i)) {
        const int lift_index = ocp.discrete().liftIndexAfterTimeStage(i);
        const bool sto = ocp.discrete().isSTOEnabledLift(lift_index);
        RiccatiFactorizer::computeCostateDirection(factorization[i], d[i], sto);
        ocp[i].expandPrimal(ocp.discrete().dt(i), s[i], d[i], sto);
      }
      else if (ocp.discrete().isTimeStageBeforeImpulse(i+1)) {
        const int impulse_index = ocp.discrete().impulseIndexAfterTimeStage(i+1);
        const bool sto = ocp.discrete().isSTOEnabledImpulse(impulse_index);
        RiccatiFactorizer::computeCostateDirection(factorization[i], d[i], sto);
        ocp[i].expandPrimal(ocp.discrete().dt(i), s[i], d[i], sto);
        d[i].setImpulseStatusByDimension(s[i].dimi());
        RiccatiFactorizer::computeLagrangeMultiplierDirection(
            factorization.switching[impulse_index], d[i], sto);
      }
      else {
        constexpr bool sto = false;
        RiccatiFactorizer::computeCostateDirection(factorization[i], d[i], sto);
        ocp[i].expandPrimal(ocp.discrete().dt(i), s[i], d[i], sto);
      }
      max_primal_step_sizes_.coeffRef(i) = ocp[i].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(i) = ocp[i].maxDualStepSize();
    }
    else if (i == N) {
      constexpr bool sto = false;
      RiccatiFactorizer::computeCostateDirection(factorization[N], d[N], sto);
      max_primal_step_sizes_.coeffRef(N) = ocp.terminal.maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(N) = ocp.terminal.maxDualStepSize();
    }
    else if (i < N + 1 + N_impulse) {
      const int impulse_index = i - (N+1);
      const bool sto = ocp.discrete().isSTOEnabledImpulse(impulse_index);
      RiccatiFactorizer::computeCostateDirection(factorization.impulse[impulse_index], 
                                                 d.impulse[impulse_index], sto);
      ocp.impulse[impulse_index].expandPrimal(s.impulse[impulse_index], 
                                              d.impulse[impulse_index]);
      max_primal_step_sizes_.coeffRef(i) 
          = ocp.impulse[impulse_index].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(i) 
          = ocp.impulse[impulse_index].maxDualStepSize();
    }
    else if (i < N + 1 + 2*N_impulse) {
      const int impulse_index = i - (N+1+N_impulse);
      const bool sto = ocp.discrete().isSTOEnabledImpulse(impulse_index);
      RiccatiFactorizer::computeCostateDirection(factorization.aux[impulse_index], 
                                                 d.aux[impulse_index], sto);
      ocp.aux[impulse_index].expandPrimal(ocp.discrete().dt_aux(impulse_index),
                                          s.aux[impulse_index], 
                                          d.aux[impulse_index], sto);
      max_primal_step_sizes_.coeffRef(i) 
          = ocp.aux[impulse_index].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(i) 
          = ocp.aux[impulse_index].maxDualStepSize();
    }
    else {
      const int lift_index = i - (N+1+2*N_impulse);
      const bool sto = ocp.discrete().isSTOEnabledLift(lift_index);
      RiccatiFactorizer::computeCostateDirection(factorization.lift[lift_index], 
                                                 d.lift[lift_index], sto);
      ocp.lift[lift_index].expandPrimal(ocp.discrete().dt_lift(lift_index),
                                        s.lift[lift_index], 
                                        d.lift[lift_index], sto);
      max_primal_step_sizes_.coeffRef(i) 
          = ocp.lift[lift_index].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(i) 
          = ocp.lift[lift_index].maxDualStepSize();
    }
  }
  N_all_ = N_all;
}


double RiccatiRecursion::maxPrimalStepSize() const {
  return max_primal_step_sizes_.head(N_all_).minCoeff();
}


double RiccatiRecursion::maxDualStepSize() const {
  return max_dual_step_sizes_.head(N_all_).minCoeff();
}


void RiccatiRecursion::getStateFeedbackGain(const int time_stage, 
                                            Eigen::MatrixXd& Kq, 
                                            Eigen::MatrixXd& Kv) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  Kq = lqr_policy_[time_stage].Kq();
  Kv = lqr_policy_[time_stage].Kv();
}

} // namespace robotoc