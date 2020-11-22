#ifndef IDOCP_DIRECTION_CALCULATOR_HXX_ 
#define IDOCP_DIRECTION_CALCULATOR_HXX_

#include "idocp/ocp/ocp_direction_calculator.hpp"
#include "idocp/ocp/riccati_factorizer.hpp"
#include "idocp/impulse/impulse_riccati_factorizer.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>


namespace idocp {

inline OCPDirectionCalculator::OCPDirectionCalculator(const double T, 
                                                      const int N, 
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


inline OCPDirectionCalculator::OCPDirectionCalculator()
  : T_(0),
    dtau_(0),
    N_(0),
    num_proc_(0),
    max_primal_step_sizes_(), 
    max_dual_step_sizes_() {
}


inline OCPDirectionCalculator::~OCPDirectionCalculator() {
}


inline void OCPDirectionCalculator::computeInitialStateDirection(
    const std::vector<Robot>& robots, const Eigen::VectorXd& q, 
    const Eigen::VectorXd& v, const HybridSolution& s, HybridDirection& d) {
  assert(q.size() == robots[0].dimq());
  assert(v.size() == robots[0].dimv());
  robots[0].subtractConfiguration(q, s[0].q, d[0].dq());
  d[0].dv() = v - s[0].v;
}


inline void OCPDirectionCalculator::computeDirection(
    HybridOCP& split_ocps, std::vector<Robot>& robots, 
    const ContactSequence& contact_sequence, 
    const HybridRiccatiFactorizer& factorizer, 
    HybridRiccatiFactorization& factorization, 
    const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization, 
    const HybridSolution& s, HybridDirection& d) {
  assert(robots.size() == num_proc_);
  const int N_impulse = contact_sequence.totalNumImpulseStages();
  const int N_lift = contact_sequence.totalNumLiftStages();
  const int N_all = N_ + 1 + 2 * N_impulse + N_lift;
  const bool exist_state_constraint = contact_sequence.existImpulseStage();
  const Eigen::VectorXd& dx0 = d[0].dx();
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_all; ++i) {
    if (i == 0) {
      aggregateLagrangeMultiplierDirection(
          contact_sequence, constraint_factorization, d.impulse, 0, 
          factorization[0]);
      computePrimalDirectionInitial(factorizer[0], factorization[0], d[0], 
                                    exist_state_constraint);
      split_ocps[0].computeCondensedPrimalDirection(robots[omp_get_thread_num()], 
                                                    dtau(contact_sequence, 0), 
                                                    s[0], d[0]);
      max_primal_step_sizes_.coeffRef(0) = split_ocps[0].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(0) = split_ocps[0].maxDualStepSize();
    }
    else if (i < N_) {
      aggregateLagrangeMultiplierDirection(
          contact_sequence, constraint_factorization, d.impulse, i, 
          factorization[i]);
      computePrimalDirection(factorizer[i], factorization[i], dx0, d[i], 
                             exist_state_constraint);
      split_ocps[i].computeCondensedPrimalDirection(robots[omp_get_thread_num()], 
                                                    dtau(contact_sequence, i), 
                                                    s[i], d[i]);
      max_primal_step_sizes_.coeffRef(i) = split_ocps[i].maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(i) = split_ocps[i].maxDualStepSize();
    }
    else if (i == N_) {
      computePrimalDirectionTerminal(factorization[N_], dx0, d[N_]);
      max_primal_step_sizes_.coeffRef(N_) = split_ocps.terminal.maxPrimalStepSize();
      max_dual_step_sizes_.coeffRef(N_) = split_ocps.terminal.maxDualStepSize();
    }
    else if (i < N_ + 1 + N_impulse) {
      const int impulse_index  = i - (N_+1);
      aggregateLagrangeMultiplierDirectionImpulse(
          contact_sequence, constraint_factorization, d.impulse, impulse_index, 
          factorization.impulse[impulse_index]);
      computePrimalDirectionImpulse(factorization.impulse[impulse_index], 
                                    dx0, d.impulse[impulse_index]);
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
      aggregateLagrangeMultiplierDirectionAux(
          contact_sequence, constraint_factorization, d.impulse, impulse_index, 
          factorization.aux[impulse_index]);
      computePrimalDirection(factorizer.aux[impulse_index], 
                             factorization.aux[impulse_index], dx0, 
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
      aggregateLagrangeMultiplierDirectionLift(
          contact_sequence, constraint_factorization, d.impulse, lift_index, 
          factorization.lift[lift_index]);
      computePrimalDirection(factorizer.lift[lift_index], 
                             factorization.lift[lift_index], dx0, 
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


inline double OCPDirectionCalculator::maxPrimalStepSize(
    const ContactSequence& contact_sequence) const {
  const int N_impulse = contact_sequence.totalNumImpulseStages();
  const int N_lift = contact_sequence.totalNumLiftStages();
  const int N_all = N_ + 1 + 2 * N_impulse + N_lift;
  return max_primal_step_sizes_.head(N_all).minCoeff();
}


inline double OCPDirectionCalculator::maxDualStepSize(
    const ContactSequence& contact_sequence) const {
  const int N_impulse = contact_sequence.totalNumImpulseStages();
  const int N_lift = contact_sequence.totalNumLiftStages();
  const int N_all = N_ + 1 + 2 * N_impulse + N_lift;
  return max_dual_step_sizes_.head(N_all).minCoeff();
}


inline void OCPDirectionCalculator::computePrimalDirectionInitial(
    const RiccatiFactorizer factorizer, 
    const RiccatiFactorization factorization, SplitDirection& d, 
    const bool exist_state_constraint) {
  RiccatiFactorizer::computeCostateDirection(factorization, d, 
                                             exist_state_constraint);
  factorizer.computeControlInputDirection(factorization, d, 
                                          exist_state_constraint);
}


template <typename VectorType>
inline void OCPDirectionCalculator::computePrimalDirection(
    const RiccatiFactorizer factorizer, 
    const RiccatiFactorization factorization, 
    const Eigen::MatrixBase<VectorType>& dx0, SplitDirection& d, 
    const bool exist_state_constraint) {
  RiccatiFactorizer::computeStateDirection(factorization, dx0, d,
                                           exist_state_constraint);
  RiccatiFactorizer::computeCostateDirection(factorization, d, 
                                             exist_state_constraint);
  factorizer.computeControlInputDirection(factorization, d, 
                                          exist_state_constraint);
}


template <typename VectorType>
inline void OCPDirectionCalculator::computePrimalDirectionTerminal(
    const RiccatiFactorization factorization, 
    const Eigen::MatrixBase<VectorType>& dx0, SplitDirection& d) {
  RiccatiFactorizer::computeStateDirection(factorization, dx0, d, false);
  RiccatiFactorizer::computeCostateDirection(factorization, d, false);
}
 

template <typename VectorType>
inline void OCPDirectionCalculator::computePrimalDirectionImpulse(
    const RiccatiFactorization factorization, 
    const Eigen::MatrixBase<VectorType>& dx0, ImpulseSplitDirection& d) {
  ImpulseRiccatiFactorizer::computeStateDirection(factorization, dx0, d);
  ImpulseRiccatiFactorizer::computeCostateDirection(factorization, d);
}
 

inline void OCPDirectionCalculator::aggregateLagrangeMultiplierDirection(
    const ContactSequence& contact_sequence,
    const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
    const std::vector<ImpulseSplitDirection>& d_impulse, const int time_stage,
    RiccatiFactorization& riccati_factorization) {
  assert(time_stage >= 0);
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  riccati_factorization.n.setZero();
  for (int i=num_impulse-1; i>=0; --i) {
    if (contact_sequence.timeStageBeforeImpulse(i) < time_stage) {
      break;
    }
    else {
      riccati_factorization.n.noalias() 
          += constraint_factorization[i].T(time_stage) * d_impulse[i].dxi();
    }
  }
}


inline void OCPDirectionCalculator::aggregateLagrangeMultiplierDirectionImpulse(
    const ContactSequence& contact_sequence,
    const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
    const std::vector<ImpulseSplitDirection>& d_impulse, 
    const int impulse_index,
    RiccatiFactorization& impulse_riccati_factorization) {
  assert(impulse_index >= 0);
  assert(impulse_index < contact_sequence.totalNumImpulseStages());
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  impulse_riccati_factorization.n.setZero();
  for (int i=num_impulse-1; i>=impulse_index; --i) {
    impulse_riccati_factorization.n.noalias() 
        += constraint_factorization[i].T_impulse(impulse_index) * d_impulse[i].dxi();
  }
}


inline void OCPDirectionCalculator::aggregateLagrangeMultiplierDirectionAux(
    const ContactSequence& contact_sequence,
    const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
    const std::vector<ImpulseSplitDirection>& d_impulse, 
    const int impulse_index,
    RiccatiFactorization& aux_riccati_factorization) {
  assert(impulse_index >= 0);
  assert(impulse_index < contact_sequence.totalNumImpulseStages());
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  aux_riccati_factorization.n.setZero();
  for (int i=num_impulse-1; i>=impulse_index; --i) {
    aux_riccati_factorization.n.noalias() 
        += constraint_factorization[i].T_aux(impulse_index) * d_impulse[i].dxi();
  }
}


inline void OCPDirectionCalculator::aggregateLagrangeMultiplierDirectionLift(
    const ContactSequence& contact_sequence,
    const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
    const std::vector<ImpulseSplitDirection>& d_impulse, 
    const int lift_index,
    RiccatiFactorization& lift_riccati_factorization) {
  assert(lift_index >= 0);
  assert(lift_index < contact_sequence.totalNumLiftStages());
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  const int time_stage_before_lift = contact_sequence.timeStageBeforeLift(lift_index);
  lift_riccati_factorization.n.setZero();
  for (int i=num_impulse-1; i>=0; --i) {
    if (contact_sequence.timeStageBeforeImpulse(i) < time_stage_before_lift) {
      break;
    }
    else {
      lift_riccati_factorization.n.noalias() 
          += constraint_factorization[i].T_lift(lift_index) * d_impulse[i].dxi();
    }
  }
}


inline double OCPDirectionCalculator::dtau(
    const ContactSequence& contact_sequence, const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  if (contact_sequence.existImpulseStage(time_stage)) {
    return (contact_sequence.impulseTime(contact_sequence.impulseIndex(time_stage)) 
              - time_stage * dtau_);
  }
  else if (contact_sequence.existLiftStage(time_stage)) {
    return (contact_sequence.liftTime(contact_sequence.liftIndex(time_stage)) 
              - time_stage * dtau_);
  }
  else {
    return dtau_;
  }
}

} // namespace idocp 

#endif // IDOCP_OCP_LINEARIZER_HXX_ 