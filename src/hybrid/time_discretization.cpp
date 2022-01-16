#include "robotoc/hybrid/time_discretization.hpp"

#include <iomanip>


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
  os << " -----------------------------------------------------------------------" << std::endl;
  os << "  grid point | grid count |      t |     dt | phase |  sto  | sto_next |" << std::endl;
  os << " -----------------------------------------------------------------------" << std::endl;
  for (int i=0; i<N(); ++i) {
    os << "  stage: " << std::right << std::setw(3) << i << " | "
       << "       " << std::right << std::setw(3) << timeStageInPhase(i) << " | " 
       << std::fixed << std::setprecision(4) << t(i) << " | " << dt(i)
       << " |   " << std::setw(3) << contactPhase(i)
       << " | " << std::setw(5) << std::boolalpha << isSTOEnabledPhase(contactPhase(i)) 
       << " |   " << std::setw(5) << std::boolalpha << isSTOEnabledNextPhase(contactPhase(i)) 
       << "  |" << std::endl;
    if (isTimeStageBeforeImpulse(i)) {
      const int impulse_index = impulseIndexAfterTimeStage(i);
      os << "    aux: " << std::right << std::setw(3) << impulse_index << " | "
         << std::right << std::setw(3) << "         0" << " | " 
         << std::fixed << std::setprecision(4) << t_impulse(impulse_index) 
         << " | " << dt_aux(impulse_index)
         << " |   " << std::setw(3) << contactPhaseAfterImpulse(impulse_index)
         << " | " << std::setw(5) << std::boolalpha 
         << isSTOEnabledPhase(contactPhaseAfterImpulse(impulse_index)) 
         << " |   " << std::setw(5) << std::boolalpha 
         << isSTOEnabledNextPhase(contactPhaseAfterImpulse(impulse_index)) 
         << "  |" << std::endl;
    }
    else if (isTimeStageBeforeLift(i)) {
      const int lift_index = liftIndexAfterTimeStage(i);
      os << "   lift: " << std::right << std::setw(3) << lift_index << " | "
         << std::right << std::setw(3) << "         0" << " | " 
         << std::fixed << std::setprecision(4) << t_lift(lift_index) 
         << " | " << dt_lift(lift_index) 
         << " |   " << std::setw(3) << contactPhaseAfterLift(lift_index)
         << " | " << std::setw(5) << std::boolalpha
         << isSTOEnabledPhase(contactPhaseAfterLift(lift_index)) 
         << " |   " << std::setw(5) << std::boolalpha
         << isSTOEnabledNextPhase(contactPhaseAfterLift(lift_index))
         << "  |" << std::endl;
    }
  }
  os << "  stage: " << std::right << std::setw(3) << N() << " | " 
     << "       " << std::right << std::setw(3) << timeStageInPhase(N()) << " | " 
     << std::fixed << std::setprecision(4) << t(N()) 
     << " |        |       |       |          |" << std::endl; 
  os << " -----------------------------------------------------------------------" << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const TimeDiscretization& discretization) {
  discretization.disp(os);
  return os;
}

} // namespace robotoc 