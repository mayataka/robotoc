#ifndef IDOCP_RICCATI_DIRECTION_CALCULATOR_HXX_ 
#define IDOCP_RICCATI_DIRECTION_CALCULATOR_HXX_

#include "idocp/ocp/riccati_direction_calculator.hpp"

#include <cassert>

namespace idocp {

inline double RiccatiDirectionCalculator::dtau(
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


inline const SplitRiccatiFactorization& RiccatiDirectionCalculator::
next_riccati_factorization(const RiccatiFactorization& factorization, 
                           const ContactSequence& contact_sequence, 
                           const int time_stage) {
  if (contact_sequence.existImpulseStage(time_stage)) {
    return factorization.impulse[contact_sequence.impulseIndex(time_stage)];
  }
  else if (contact_sequence.existLiftStage(time_stage)) {
    return factorization.lift[contact_sequence.liftIndex(time_stage)];
  }
  else {
    return factorization[time_stage+1];
  }
}

} // namespace idocp 

#endif // IDOCP_RICCATI_DIRECTION_CALCULATOR_HXX_ 