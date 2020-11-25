#ifndef IDOCP_DIRECTION_CALCULATOR_HXX_ 
#define IDOCP_DIRECTION_CALCULATOR_HXX_

#include "idocp/ocp/ocp_direction_calculator.hpp"
#include "idocp/ocp/riccati_factorizer.hpp"
#include "idocp/impulse/impulse_riccati_factorizer.hpp"

#include <cassert>

namespace idocp {

inline void OCPDirectionCalculator::computePrimalDirectionInitial(
    const RiccatiFactorizer factorizer, 
    const RiccatiFactorization factorization, 
    const RiccatiFactorization factorization_next, SplitDirection& d, 
    const bool exist_state_constraint) {
  RiccatiFactorizer::computeCostateDirection(factorization, d, 
                                             exist_state_constraint);
  factorizer.computeControlInputDirection(factorization_next, d, 
                                          exist_state_constraint);
}


template <typename VectorType>
inline void OCPDirectionCalculator::computePrimalDirection(
    const RiccatiFactorizer factorizer, 
    const RiccatiFactorization factorization, 
    const RiccatiFactorization factorization_next, 
    const Eigen::MatrixBase<VectorType>& dx0, SplitDirection& d, 
    const bool exist_state_constraint) {
  RiccatiFactorizer::computeStateDirection(factorization, dx0, d,
                                           exist_state_constraint);
  RiccatiFactorizer::computeCostateDirection(factorization, d, 
                                             exist_state_constraint);
  factorizer.computeControlInputDirection(factorization_next, d, 
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
    const Eigen::MatrixBase<VectorType>& dx0, ImpulseSplitDirection& d,
    const bool exist_state_constraint) {
  ImpulseRiccatiFactorizer::computeStateDirection(factorization, dx0, d, 
                                                  exist_state_constraint);
  ImpulseRiccatiFactorizer::computeCostateDirection(factorization, d, 
                                                    exist_state_constraint);
}
 

inline void OCPDirectionCalculator::aggregateLagrangeMultiplierDirection(
    const ContactSequence& contact_sequence,
    const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
    const std::vector<ImpulseSplitDirection>& d_impulse, const int time_stage,
    RiccatiFactorization& riccati_factorization) const {
  assert(time_stage >= 0);
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  for (int i=time_stage; i<N_; ++i) {
    if (contact_sequence.existImpulseStage(i)) {
      const int impulse_index = contact_sequence.impulseIndex(i);
      riccati_factorization.n.noalias() += 
        constraint_factorization[impulse_index].T(time_stage) * d_impulse[impulse_index].dxi();
    }
  }
}


inline void OCPDirectionCalculator::aggregateLagrangeMultiplierDirectionImpulse(
    const ContactSequence& contact_sequence,
    const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
    const std::vector<ImpulseSplitDirection>& d_impulse, 
    const int impulse_index,
    RiccatiFactorization& impulse_riccati_factorization) const {
  assert(impulse_index >= 0);
  assert(impulse_index < contact_sequence.totalNumImpulseStages());
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  impulse_riccati_factorization.n.setZero();
  for (int i=impulse_index; i<num_impulse; ++i) {
    impulse_riccati_factorization.n.noalias() 
        += constraint_factorization[i].T_impulse(impulse_index) * d_impulse[i].dxi();
  }
}


inline void OCPDirectionCalculator::aggregateLagrangeMultiplierDirectionAux(
    const ContactSequence& contact_sequence,
    const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
    const std::vector<ImpulseSplitDirection>& d_impulse, 
    const int impulse_index,
    RiccatiFactorization& aux_riccati_factorization) const {
  assert(impulse_index >= 0);
  assert(impulse_index < contact_sequence.totalNumImpulseStages());
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  aux_riccati_factorization.n.setZero();
  for (int i=impulse_index+1; i<num_impulse; ++i) {
    aux_riccati_factorization.n.noalias() 
        += constraint_factorization[i].T_aux(impulse_index) * d_impulse[i].dxi();
  }
}


inline void OCPDirectionCalculator::aggregateLagrangeMultiplierDirectionLift(
    const ContactSequence& contact_sequence,
    const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
    const std::vector<ImpulseSplitDirection>& d_impulse, 
    const int lift_index,
    RiccatiFactorization& lift_riccati_factorization) const {
  assert(lift_index >= 0);
  assert(lift_index < contact_sequence.totalNumLiftStages());
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  const int time_stage_after_lift = contact_sequence.timeStageAfterLift(lift_index);
  lift_riccati_factorization.n.setZero();
  for (int i=time_stage_after_lift; i<N_; ++i) {
    if (contact_sequence.existImpulseStage(i)) {
      const int impulse_index = contact_sequence.impulseIndex(i);
      lift_riccati_factorization.n.noalias() += 
        constraint_factorization[impulse_index].T_lift(lift_index) * d_impulse[impulse_index].dxi();
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