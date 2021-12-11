#include "robotoc/hybrid/time_discretization.hpp"


namespace robotoc {

void TimeDiscretization::disp(std::ostream& os) const {
  os << "Time discretization of optimal control problem (OCP):" << std::endl;
  os << "  T: " << T_ << std::endl;
  os << "  N_ideal: " << N_ideal() << std::endl;
  os << "  N: " << N() << std::endl;
  os << "  N_impulse: " << N_impulse() << std::endl;
  os << "  N_lift: " << N_lift() << std::endl;
  os << "  N_all: " << (N()+1+2*N_impulse()+N_lift()) << std::endl;
  os << "  No. of discrete events: " << numDiscreteEvents() << std::endl;
  os << "  No. of contact phases: " << numContactPhases() << std::endl;
  for (int i=0; i<numContactPhases(); ++i) {
    os << "    No. of grids at contact phase " << i << ": " << N_phase(i) << std::endl;
  }
  os << "  isFormulationTractable: " << std::boolalpha 
     << isFormulationTractable() << std::endl;;
  os << "  isSwitchingTimeConsistent: " << std::boolalpha
     << isSwitchingTimeConsistent() << std::endl;
  for (int i=0; i<N(); ++i) {
    os << "  i = " << i << ": t = " << t(i) << ": dt = " << dt(i) 
       << ": phase = " << contactPhase(i) 
       << ": sto = " << std::boolalpha << isSTOEnabledPhase(contactPhase(i)) 
       << ": sto_next = " << std::boolalpha << isSTOEnabledNextPhase(contactPhase(i)) << std::endl;
    if (isTimeStageBeforeImpulse(i)) {
      const int impulse_index = impulseIndexAfterTimeStage(i);
      os << "  aux = " << impulse_index
         << ": t =  " << t_impulse(impulse_index) 
         << ": dt_aux = " << dt_aux(impulse_index)
         << ": phase = " << contactPhaseAfterImpulse(impulse_index)
         << ": sto = " << std::boolalpha 
         << isSTOEnabledPhase(contactPhaseAfterImpulse(impulse_index)) 
         << ": sto_next = " << std::boolalpha 
         << isSTOEnabledNextPhase(contactPhaseAfterImpulse(impulse_index)) 
         << std::endl;
    }
    else if (isTimeStageBeforeLift(i)) {
      const int lift_index = liftIndexAfterTimeStage(i);
      os << "  lift = " << lift_index 
         << ": t =  " << t_lift(lift_index) 
         << ": dt_lift = " << dt_lift(lift_index) 
         << ": phase = " << contactPhaseAfterLift(lift_index)
         << ": sto = " << std::boolalpha
         << isSTOEnabledPhase(contactPhaseAfterLift(lift_index)) 
         << ": sto_next = " << std::boolalpha
         << isSTOEnabledNextPhase(contactPhaseAfterLift(lift_index)) << std::endl;
    }
  }
  os << "  i = " << N() << ": t = " << t(N()) << std::flush; 
}


std::ostream& operator<<(std::ostream& os, 
                         const TimeDiscretization& discretization) {
  discretization.disp(os);
  return os;
}

} // namespace robotoc 