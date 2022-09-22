#include "robotoc/riccati/riccati_recursion.hpp"

#include <omp.h>
#include <stdexcept>
#include <iostream>
#include <cassert>

namespace robotoc {

RiccatiRecursion::RiccatiRecursion(const OCP& ocp, const int nthreads, 
                                   const double max_dts0)
  : nthreads_(nthreads),
    N_all_(ocp.N()+1),
    factorizer_(ocp.robot(), max_dts0),
    lqr_policy_(ocp.robot(), ocp.N(), ocp.reservedNumDiscreteEvents()),
    sto_policy_(2*ocp.reservedNumDiscreteEvents()+1, 
                STOPolicy(ocp.robot())),
    factorization_m_(ocp.robot()),
    max_primal_step_sizes_(
        Eigen::VectorXd::Zero(ocp.N()+1+3*ocp.reservedNumDiscreteEvents())), 
    max_dual_step_sizes_(
        Eigen::VectorXd::Zero(ocp.N()+1+3*ocp.reservedNumDiscreteEvents())) {
  if (nthreads <= 0) {
    throw std::out_of_range("[RiccatiRecursion] invalid argument: 'nthreads' must be positive!");
  }
}


RiccatiRecursion::RiccatiRecursion()
  : nthreads_(0),
    N_all_(0),
    factorizer_(),
    lqr_policy_(),
    sto_policy_(),
    factorization_m_(),
    max_primal_step_sizes_(), 
    max_dual_step_sizes_() {
}


RiccatiRecursion::~RiccatiRecursion() {
}


void RiccatiRecursion::setRegularization(const double max_dts0) {
  assert(max_dts0 > 0);
  factorizer_.setRegularization(max_dts0);
}


void RiccatiRecursion::reserve(const OCP& ocp) {
  const int reserved_num_discrete_events = ocp.reservedNumDiscreteEvents();
  lqr_policy_.reserve(ocp.robot(), reserved_num_discrete_events);
  while (sto_policy_.size() < 2*reserved_num_discrete_events+1) {
    sto_policy_.emplace_back(ocp.robot());
  }
  const int N = ocp.timeDiscretization().N();
  const int max_N_all = N + 1 + 3*reserved_num_discrete_events;
  if (max_N_all > max_primal_step_sizes_.size()) {
    max_primal_step_sizes_.resize(max_N_all);
    max_dual_step_sizes_.resize(max_N_all);
  }
}


void RiccatiRecursion::backwardRiccatiRecursion(
    const OCP& ocp, KKTMatrix& kkt_matrix, KKTResidual& kkt_residual, 
    RiccatiFactorization& factorization) {
  const int N = ocp.timeDiscretization().N();
  factorization[N].P = kkt_matrix[N].Qxx;
  factorization[N].s = - kkt_residual[N].lx;
  for (int i=N-1; i>=0; --i) {
    if (ocp.timeDiscretization().isTimeStageBeforeImpulse(i)) {
      assert(!ocp.timeDiscretization().isTimeStageBeforeImpulse(i+1));
      const int impulse_index = ocp.timeDiscretization().impulseIndexAfterTimeStage(i);
      const int phase = ocp.timeDiscretization().contactPhase(i);
      const bool sto = ocp.timeDiscretization().isSTOEnabledPhase(phase);
      const bool sto_next = ocp.timeDiscretization().isSTOEnabledNextPhase(phase);
      const bool sto_next_next = ocp.timeDiscretization().isSTOEnabledNextPhase(phase+1);
      factorizer_.backwardRiccatiRecursion(factorization[i+1], 
                                           kkt_matrix.aux[impulse_index], 
                                           kkt_residual.aux[impulse_index], 
                                           factorization.aux[impulse_index], 
                                           lqr_policy_.aux[impulse_index],
                                           sto_next, sto_next_next);
      if (sto || sto_next) {
        const int event_index = ocp.timeDiscretization().eventIndexImpulse(impulse_index);
        factorizer_.backwardRiccatiRecursionPhaseTransition(
            factorization.aux[impulse_index], factorization_m_,
            sto_policy_[event_index+1], sto_next_next);
        factorizer_.backwardRiccatiRecursion(factorization_m_, 
                                             kkt_matrix.impulse[impulse_index], 
                                             kkt_residual.impulse[impulse_index], 
                                             factorization.impulse[impulse_index],
                                             sto_next);
      }
      else {
        factorizer_.backwardRiccatiRecursion(factorization.aux[impulse_index], 
                                             kkt_matrix.impulse[impulse_index], 
                                             kkt_residual.impulse[impulse_index], 
                                             factorization.impulse[impulse_index],
                                             sto_next);
      }
      factorizer_.backwardRiccatiRecursion(factorization.impulse[impulse_index], 
                                           kkt_matrix[i], kkt_residual[i], 
                                           factorization[i], lqr_policy_[i],
                                           sto, sto_next);
      if (i-1 >= 0) {
        factorizer_.backwardRiccatiRecursion(factorization[i], kkt_matrix[i-1], 
                                             kkt_residual[i-1],
                                             kkt_matrix.switching[impulse_index], 
                                             kkt_residual.switching[impulse_index], 
                                             factorization[i-1], 
                                             factorization.switching[impulse_index], 
                                             lqr_policy_[i-1], sto, sto_next);
      }
    }
    else if (ocp.timeDiscretization().isTimeStageBeforeLift(i)) {
      assert(!ocp.timeDiscretization().isTimeStageBeforeImpulse(i+1));
      const int lift_index = ocp.timeDiscretization().liftIndexAfterTimeStage(i);
      const int phase = ocp.timeDiscretization().contactPhase(i);
      const bool sto = ocp.timeDiscretization().isSTOEnabledPhase(phase);
      const bool sto_next = ocp.timeDiscretization().isSTOEnabledNextPhase(phase);
      const bool sto_next_next = ocp.timeDiscretization().isSTOEnabledNextPhase(phase+1);
      factorizer_.backwardRiccatiRecursion(factorization[i+1], 
                                           kkt_matrix.lift[lift_index], 
                                           kkt_residual.lift[lift_index], 
                                           factorization.lift[lift_index], 
                                           lqr_policy_.lift[lift_index],
                                           sto_next, sto_next_next);
      if (sto || sto_next) {
        const int event_index = ocp.timeDiscretization().eventIndexLift(lift_index);
        factorizer_.backwardRiccatiRecursionPhaseTransition(
            factorization.lift[lift_index], factorization_m_,
            sto_policy_[event_index+1], sto_next_next);
        factorizer_.backwardRiccatiRecursion(factorization_m_, kkt_matrix[i], 
                                             kkt_residual[i], factorization[i], 
                                             lqr_policy_[i], sto, sto_next);
      }
      else {
        factorizer_.backwardRiccatiRecursion(factorization.lift[lift_index], 
                                             kkt_matrix[i], kkt_residual[i], 
                                             factorization[i], lqr_policy_[i], 
                                             sto, sto_next);
      }
    }
    else if (!ocp.timeDiscretization().isTimeStageBeforeImpulse(i+1)) {
      const int phase = ocp.timeDiscretization().contactPhase(i);
      const bool sto = ocp.timeDiscretization().isSTOEnabledPhase(phase);
      const bool sto_next = ocp.timeDiscretization().isSTOEnabledNextPhase(phase);
      factorizer_.backwardRiccatiRecursion(factorization[i+1], 
                                           kkt_matrix[i], kkt_residual[i], 
                                           factorization[i], lqr_policy_[i],
                                           sto, sto_next);
    }
  }
  const bool sto = ocp.timeDiscretization().isSTOEnabledPhase(0);
  if (sto) {
    const bool sto_next = ocp.timeDiscretization().isSTOEnabledNextPhase(0);
    factorizer_.backwardRiccatiRecursionPhaseTransition(factorization[0], 
                                                        factorization_m_,
                                                        sto_policy_[0], sto_next);
  }
}


void RiccatiRecursion::forwardRiccatiRecursion(
    const OCP& ocp, const KKTMatrix& kkt_matrix, 
    const KKTResidual& kkt_residual, Direction& d) const {
  d[0].dts = 0.0;
  d[0].dts_next = 0.0;
  if (ocp.timeDiscretization().isSTOEnabledPhase(0)) {
    factorizer_.computeSwitchingTimeDirection(sto_policy_[0], d[0], false);
  }
  const int N = ocp.timeDiscretization().N();
  for (int i=0; i<N; ++i) {
    if (ocp.timeDiscretization().isTimeStageBeforeImpulse(i)) {
      assert(!ocp.timeDiscretization().isTimeStageBeforeImpulse(i+1));
      const int impulse_index = ocp.timeDiscretization().impulseIndexAfterTimeStage(i);
      const int phase = ocp.timeDiscretization().contactPhase(i);
      const bool sto = ocp.timeDiscretization().isSTOEnabledPhase(phase);
      const bool sto_next = ocp.timeDiscretization().isSTOEnabledNextPhase(phase);
      const bool sto_next_next = ocp.timeDiscretization().isSTOEnabledNextPhase(phase+1);
      if (i-1 >= 0) {
        factorizer_.forwardRiccatiRecursion(kkt_matrix[i-1], kkt_residual[i-1],
                                            lqr_policy_[i-1], d[i-1], d[i], 
                                            sto, sto_next);
      }
      factorizer_.forwardRiccatiRecursion(kkt_matrix[i], kkt_residual[i],  
                                          lqr_policy_[i], d[i], 
                                          d.impulse[impulse_index],
                                          sto, sto_next);
      factorizer_.forwardRiccatiRecursion(kkt_matrix.impulse[impulse_index], 
                                          kkt_residual.impulse[impulse_index],
                                          d.impulse[impulse_index], 
                                          d.aux[impulse_index]);
      d.aux[impulse_index].dts = d.impulse[impulse_index].dts_next;
      if (sto_next_next) {
        const int next_event_index = ocp.timeDiscretization().eventIndexImpulse(impulse_index) + 1;
        factorizer_.computeSwitchingTimeDirection(sto_policy_[next_event_index], 
                                                  d.aux[impulse_index], sto_next);
      }
      else {
        d.aux[impulse_index].dts_next = 0.0;
      }
      factorizer_.forwardRiccatiRecursion(kkt_matrix.aux[impulse_index], 
                                          kkt_residual.aux[impulse_index],
                                          lqr_policy_.aux[impulse_index], 
                                          d.aux[impulse_index], d[i+1], 
                                          sto_next, sto_next_next);
    }
    else if (ocp.timeDiscretization().isTimeStageBeforeLift(i)) {
      assert(!ocp.timeDiscretization().isTimeStageBeforeImpulse(i+1));
      const int lift_index = ocp.timeDiscretization().liftIndexAfterTimeStage(i);
      const int phase = ocp.timeDiscretization().contactPhase(i);
      const bool sto = ocp.timeDiscretization().isSTOEnabledPhase(phase);
      const bool sto_next = ocp.timeDiscretization().isSTOEnabledNextPhase(phase);
      const bool sto_next_next = ocp.timeDiscretization().isSTOEnabledNextPhase(phase+1);
      factorizer_.forwardRiccatiRecursion(kkt_matrix[i], kkt_residual[i],  
                                          lqr_policy_[i], d[i], 
                                          d.lift[lift_index], sto, sto_next);
      d.lift[lift_index].dts = d[i].dts_next;
      if (sto_next_next) {
        const int next_event_index = ocp.timeDiscretization().eventIndexLift(lift_index) + 1;
        factorizer_.computeSwitchingTimeDirection(sto_policy_[next_event_index], 
                                                  d.lift[lift_index], sto_next);
      }
      else {
        d.lift[lift_index].dts_next = 0.0;
      }
      factorizer_.forwardRiccatiRecursion(kkt_matrix.lift[lift_index], 
                                          kkt_residual.lift[lift_index], 
                                          lqr_policy_.lift[lift_index], 
                                          d.lift[lift_index], d[i+1], 
                                          sto_next, sto_next_next);
    }
    else if (!ocp.timeDiscretization().isTimeStageBeforeImpulse(i+1)) {
      const int phase = ocp.timeDiscretization().contactPhase(i);
      const bool sto = ocp.timeDiscretization().isSTOEnabledPhase(phase);
      const bool sto_next = ocp.timeDiscretization().isSTOEnabledNextPhase(phase);
      factorizer_.forwardRiccatiRecursion(kkt_matrix[i], kkt_residual[i],  
                                          lqr_policy_[i], d[i], d[i+1], 
                                          sto, sto_next);
    }
  }
}


void RiccatiRecursion::computeDirection(
    OCP& ocp, const std::shared_ptr<ContactSequence>& contact_sequence, 
    const RiccatiFactorization& factorization, Direction& d) {
  const int N = ocp.timeDiscretization().N();
  const int N_impulse = ocp.timeDiscretization().N_impulse();
  const int N_lift = ocp.timeDiscretization().N_lift();
  const int N_all = N + 1 + 2*N_impulse + N_lift;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N) {
      const int phase = ocp.timeDiscretization().contactPhase(i);
      const bool sto = ocp.timeDiscretization().isSTOEnabledPhase(phase);
      const bool sto_next = ocp.timeDiscretization().isSTOEnabledNextPhase(phase);
      RiccatiFactorizer::computeCostateDirection(factorization[i], d[i], 
                                                 sto, sto_next);
      ocp[i].expandPrimal(contact_sequence->contactStatus(phase), d[i]);
      if (ocp.timeDiscretization().isTimeStageBeforeImpulse(i+1)) {
        const int impulse_index = ocp.timeDiscretization().impulseIndexAfterTimeStage(i+1);
        d[i].setSwitchingConstraintDimension(contact_sequence->impulseStatus(impulse_index).dimi());
        RiccatiFactorizer::computeLagrangeMultiplierDirection(
            factorization.switching[impulse_index], d[i], sto, sto_next);
      }
      max_primal_step_sizes_.coeffRef(i) = ocp[i].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(i) = ocp[i].maxDualStepSize();
    }
    else if (i == N) {
      constexpr bool sto = false;
      constexpr bool sto_next = false;
      RiccatiFactorizer::computeCostateDirection(factorization[N], d[N], 
                                                 sto, sto_next);
      max_primal_step_sizes_.coeffRef(N) = ocp.terminal.maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(N) = ocp.terminal.maxDualStepSize();
    }
    else if (i < N + 1 + N_impulse) {
      const int impulse_index = i - (N+1);
      const int phase = ocp.timeDiscretization().contactPhaseAfterImpulse(impulse_index);
      const bool sto = ocp.timeDiscretization().isSTOEnabledPhase(phase);
      const bool sto_next = ocp.timeDiscretization().isSTOEnabledNextPhase(phase);
      RiccatiFactorizer::computeCostateDirection(
          factorization.impulse[impulse_index], d.impulse[impulse_index], sto);
      ocp.impulse[impulse_index].expandPrimal(
          contact_sequence->impulseStatus(impulse_index), 
          d.impulse[impulse_index]);
      max_primal_step_sizes_.coeffRef(i) 
          = ocp.impulse[impulse_index].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(i) 
          = ocp.impulse[impulse_index].maxDualStepSize();
    }
    else if (i < N + 1 + 2*N_impulse) {
      const int impulse_index = i - (N+1+N_impulse);
      const int phase = ocp.timeDiscretization().contactPhaseAfterImpulse(impulse_index);
      const bool sto = ocp.timeDiscretization().isSTOEnabledPhase(phase);
      const bool sto_next = ocp.timeDiscretization().isSTOEnabledNextPhase(phase);
      RiccatiFactorizer::computeCostateDirection(factorization.aux[impulse_index], 
                                                 d.aux[impulse_index], 
                                                 sto, sto_next);
      ocp.aux[impulse_index].expandPrimal(contact_sequence->contactStatus(phase),
                                          d.aux[impulse_index]);
      max_primal_step_sizes_.coeffRef(i) 
          = ocp.aux[impulse_index].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(i) 
          = ocp.aux[impulse_index].maxDualStepSize();
    }
    else {
      const int lift_index = i - (N+1+2*N_impulse);
      const int phase = ocp.timeDiscretization().contactPhaseAfterLift(lift_index);
      const bool sto = ocp.timeDiscretization().isSTOEnabledPhase(phase);
      const bool sto_next = ocp.timeDiscretization().isSTOEnabledNextPhase(phase);
      RiccatiFactorizer::computeCostateDirection(factorization.lift[lift_index], 
                                                 d.lift[lift_index], 
                                                 sto, sto_next);
      ocp.lift[lift_index].expandPrimal(contact_sequence->contactStatus(phase),
                                        d.lift[lift_index]);
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


const hybrid_container<LQRPolicy>& RiccatiRecursion::getLQRPolicy() const {
  return lqr_policy_;
}

} // namespace robotoc