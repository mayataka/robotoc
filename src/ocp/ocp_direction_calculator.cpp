#include "idocp/ocp/ocp_direction_calculator.hpp"

#include "idocp/ocp/riccati_factorizer.hpp"
#include "idocp/impulse/impulse_riccati_factorizer.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>


namespace idocp {

OCPDirectionCalculator::OCPDirectionCalculator(const double T, const int N, 
                                               const int max_num_impulse, 
                                               const int num_proc) 
  : T_(T),
    dtau_(T/N),
    N_(N),
    num_proc_(num_proc),
    max_primal_step_sizes_(Eigen::VectorXd::Zero(N+1+3*max_num_impulse)), 
    max_dual_step_sizes_(Eigen::VectorXd::Zero(N+1+3*max_num_impulse)) {
  try {
    if (T <= 0) {
      throw std::out_of_range("invalid value: T must be positive!");
    }
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
    if (max_num_impulse < 0) {
      throw std::out_of_range("invalid value: max_num_impulse must be non-negative!");
    }
    if (num_proc <= 0) {
      throw std::out_of_range("invalid value: num_proc must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


OCPDirectionCalculator::OCPDirectionCalculator()
  : T_(0),
    dtau_(0),
    N_(0),
    num_proc_(0),
    max_primal_step_sizes_(), 
    max_dual_step_sizes_() {
}


OCPDirectionCalculator::~OCPDirectionCalculator() {
}


void OCPDirectionCalculator::computeInitialStateDirection(
    const std::vector<Robot>& robots, const Eigen::VectorXd& q, 
    const Eigen::VectorXd& v, const HybridSolution& s, HybridDirection& d) {
  assert(q.size() == robots[0].dimq());
  assert(v.size() == robots[0].dimv());
  robots[0].subtractConfiguration(q, s[0].q, d[0].dq());
  d[0].dv() = v - s[0].v;
}


// void OCPDirectionCalculator::aggregateLagrangeMultiplierDirection(
//     const ContactSequence& contact_sequence, 
//     HybridRiccatiFactorization& factorization, 
//     const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
//     const HybridDirection& d) const {
//   const int N_impulse = contact_sequence.totalNumImpulseStages();
//   const int N_lift = contact_sequence.totalNumLiftStages();
//   const int N_all = N_ + 2 * N_impulse + N_lift;
//   #pragma omp parallel for num_threads(num_proc_)
//   for (int i=0; i<N_all; ++i) {
//     if (i < N_) {
//       aggregateLagrangeMultiplierDirection(
//           contact_sequence, constraint_factorization, d.impulse, i, 
//           factorization[i]);
//     }
//     else if (i < N_ + N_impulse) {
//       const int impulse_index = i - N_;
//       aggregateLagrangeMultiplierDirectionImpulse(
//           contact_sequence, constraint_factorization, d.impulse, impulse_index, 
//           factorization.impulse[impulse_index]);
//     }
//     else if (i < N_ + 2*N_impulse) {
//       const int impulse_index = i - N_ - N_impulse;
//       aggregateLagrangeMultiplierDirectionAux(
//           contact_sequence, constraint_factorization, d.impulse, impulse_index, 
//           factorization.aux[impulse_index]);
//     }
//     else {
//       const int lift_index = i - N_ - 2*N_impulse;
//       aggregateLagrangeMultiplierDirectionLift(
//           contact_sequence, constraint_factorization, d.impulse, lift_index, 
//           factorization.lift[lift_index]);
//     }
//   }
// }


void OCPDirectionCalculator::computeDirection(
    HybridOCP& split_ocps, std::vector<Robot>& robots, 
    const ContactSequence& contact_sequence, 
    const HybridRiccatiFactorizer& factorizer, 
    HybridRiccatiFactorization& factorization, const HybridSolution& s, 
    HybridDirection& d) {
  assert(robots.size() == num_proc_);
  const int N_impulse = contact_sequence.totalNumImpulseStages();
  const int N_lift = contact_sequence.totalNumLiftStages();
  const int N_all = N_ + 1 + 2 * N_impulse + N_lift;
  const bool exist_state_constraint = contact_sequence.existImpulseStage();
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_all; ++i) {
    if (i < N_) {
      if (contact_sequence.existImpulseStage(i)) {
        computePrimalDirection(factorizer[i], factorization[i], 
                               factorization.impulse[contact_sequence.impulseIndex(i)], 
                               d[i], exist_state_constraint);
      }
      else if (contact_sequence.existLiftStage(i)) {
        computePrimalDirection(factorizer[i], factorization[i],
                               factorization.lift[contact_sequence.liftIndex(i)], 
                               d[i], exist_state_constraint);
      }
      else {
        computePrimalDirection(factorizer[i], factorization[i], 
                               factorization[i+1], d[i], exist_state_constraint);
      }
      split_ocps[i].computeCondensedPrimalDirection(robots[omp_get_thread_num()], 
                                                    dtau(contact_sequence, i), 
                                                    s[i], d[i]);
      max_primal_step_sizes_.coeffRef(i) = split_ocps[i].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(i) = split_ocps[i].maxDualStepSize();
    }
    else if (i == N_) {
      computePrimalDirectionTerminal(factorization[N_], d[N_]);
      max_primal_step_sizes_.coeffRef(N_) = split_ocps.terminal.maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(N_) = split_ocps.terminal.maxDualStepSize();
    }
    else if (i < N_ + 1 + N_impulse) {
      const int impulse_index  = i - (N_+1);
      computePrimalDirectionImpulse(factorization.impulse[impulse_index], 
                                    d.impulse[impulse_index], exist_state_constraint);
      split_ocps.impulse[impulse_index].computeCondensedPrimalDirection(
          robots[omp_get_thread_num()], s.impulse[impulse_index], 
          d.impulse[impulse_index]);
      max_primal_step_sizes_.coeffRef(i) 
          = split_ocps.impulse[impulse_index].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(i) 
          = split_ocps.impulse[impulse_index].maxDualStepSize();
    }
    else if (i < N_ + 1 + 2*N_impulse) {
      const int impulse_index  = i - (N_+1+N_impulse);
      const int time_stage_after_impulse 
          = contact_sequence.timeStageAfterImpulse(impulse_index);
      const double dtau_aux 
          = time_stage_after_impulse * dtau_ 
              - contact_sequence.impulseTime(impulse_index);
      computePrimalDirection(factorizer.aux[impulse_index], 
                             factorization.aux[impulse_index], 
                             factorization[time_stage_after_impulse], 
                             d.aux[impulse_index], exist_state_constraint);
      split_ocps.aux[impulse_index].computeCondensedPrimalDirection(
          robots[omp_get_thread_num()], dtau_aux, s.aux[impulse_index], 
          d.aux[impulse_index]);
      max_primal_step_sizes_.coeffRef(i) 
          = split_ocps.aux[impulse_index].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(i) 
          = split_ocps.aux[impulse_index].maxDualStepSize();
    }
    else {
      const int lift_index = i - (N_+1+2*N_impulse);
      const int time_stage_after_lift 
          = contact_sequence.timeStageAfterLift(lift_index);
      const double dtau_aux
          = time_stage_after_lift * dtau_ 
              - contact_sequence.liftTime(lift_index);
      computePrimalDirection(factorizer.lift[lift_index], 
                             factorization.lift[lift_index], 
                             factorization[time_stage_after_lift], 
                             d.lift[lift_index], exist_state_constraint);
      split_ocps.lift[lift_index].computeCondensedPrimalDirection(
          robots[omp_get_thread_num()], dtau_aux, s.lift[lift_index], 
          d.lift[lift_index]);
      max_primal_step_sizes_.coeffRef(i) 
          = split_ocps.lift[lift_index].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(i) 
          = split_ocps.lift[lift_index].maxDualStepSize();
    }
  }
}


double OCPDirectionCalculator::maxPrimalStepSize(
    const ContactSequence& contact_sequence) const {
  const int N_impulse = contact_sequence.totalNumImpulseStages();
  const int N_lift = contact_sequence.totalNumLiftStages();
  const int N_all = N_ + 1 + 2 * N_impulse + N_lift;
  return max_primal_step_sizes_.head(N_all).minCoeff();
}


double OCPDirectionCalculator::maxDualStepSize(
    const ContactSequence& contact_sequence) const {
  const int N_impulse = contact_sequence.totalNumImpulseStages();
  const int N_lift = contact_sequence.totalNumLiftStages();
  const int N_all = N_ + 1 + 2 * N_impulse + N_lift;
  return max_dual_step_sizes_.head(N_all).minCoeff();
}

} // namespace idocp 