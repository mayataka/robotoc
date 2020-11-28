#ifndef IDOCP_DIRECTION_CALCULATOR_HXX_ 
#define IDOCP_DIRECTION_CALCULATOR_HXX_

#include "idocp/ocp/ocp_direction_calculator.hpp"
#include "idocp/ocp/riccati_factorizer.hpp"
#include "idocp/impulse/impulse_riccati_factorizer.hpp"

#include <cassert>

namespace idocp {

inline void OCPDirectionCalculator::computePrimalDirection(
    const RiccatiFactorizer factorizer, 
    const RiccatiFactorization factorization, 
    const RiccatiFactorization factorization_next, SplitDirection& d, 
    const bool exist_state_constraint) {
  RiccatiFactorizer::computeCostateDirection(factorization, d, 
                                             exist_state_constraint);
  factorizer.computeControlInputDirection(factorization_next, d, 
                                          exist_state_constraint);
}


inline void OCPDirectionCalculator::computePrimalDirectionTerminal(
    const RiccatiFactorization factorization, SplitDirection& d) {
  RiccatiFactorizer::computeCostateDirection(factorization, d, false);
}
 

inline void OCPDirectionCalculator::computePrimalDirectionImpulse(
    const RiccatiFactorization factorization, ImpulseSplitDirection& d,
    const bool exist_state_constraint) {
  ImpulseRiccatiFactorizer::computeCostateDirection(factorization, d, 
                                                    exist_state_constraint);
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