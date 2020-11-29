#ifndef IDOCP_DIRECTION_CALCULATOR_HXX_ 
#define IDOCP_DIRECTION_CALCULATOR_HXX_

#include "idocp/ocp/ocp_direction_calculator.hpp"
#include "idocp/ocp/riccati_factorizer.hpp"
#include "idocp/impulse/impulse_riccati_factorizer.hpp"

#include <cassert>

namespace idocp {

inline const RiccatiFactorization& OCPDirectionCalculator::
riccati_factorization_next(const HybridRiccatiFactorization& factorization, 
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