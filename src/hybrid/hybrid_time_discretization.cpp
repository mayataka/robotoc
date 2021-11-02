#include "robotoc/hybrid/hybrid_ocp_discretization.hpp"


namespace robotoc {

void HybridOCPDiscretization::disp(std::ostream& os) const {
  os << "The discretized optimal control problem (OCP):" << std::endl;
  os << "  T = " << T_ << std::endl;
  os << "  N_ideal = " << N_ideal() << std::endl;
  os << "  N = " << N() << std::endl;
  os << "  N_impulse = " << N_impulse() << std::endl;
  os << "  N_lift = " << N_lift() << std::endl;
  os << "  N_all = " << (N()+1+2*N_impulse()+N_lift()) << std::endl;
  os << "  No. of discrete events = " << numDiscreteEvents() << std::endl;
  os << "  No. of contact phases = " << numContactPhases() << std::endl;
  for (int i=0; i<numContactPhases(); ++i) {
    os << "    No. of grids at contact phase: " << i << " = " << N_phase(i) << std::endl;
  }
  os << "  isFormulationTractable: ";
  if (isFormulationTractable()) os << "true" << std::endl;
  else os << "false" << std::endl;
  os << "  isSwitchingTimeConsistent: ";
  if (isSwitchingTimeConsistent()) os << "true" << std::endl;
  else os << "false" << std::endl;
  for (int i=0; i<N(); ++i) {
    os << "  i = " << i << ": t = " << t(i) << ": dt = " << dt(i) 
       << ": phase = " << contactPhase(i) << std::endl; 
    if (isTimeStageBeforeImpulse(i)) {
      const int impulse_index = impulseIndexAfterTimeStage(i);
      os << "  impulse = " << impulse_index
                << ": t =  " << t_impulse(impulse_index) 
                << ": dt_aux = " << dt_aux(impulse_index)
                << ": sto = ";
      if (isSTOEnabledImpulse(impulse_index)) os << "true" << std::endl;
      else os << "false" << std::endl;
    }
    else if (isTimeStageBeforeLift(i)) {
      const int lift_index = liftIndexAfterTimeStage(i);
      os << "  lift = " << lift_index 
                << ": t =  " << t_lift(lift_index) 
                << ": dt_lift = " << dt_lift(lift_index) << std::endl;
      if (isSTOEnabledLift(lift_index)) os << "true" << std::endl;
      else os << "false" << std::endl;
    }
  }
  os << "  i = " << N() << ": t = " << t(N()) << std::flush; 
}


std::ostream& operator<<(std::ostream& os, 
                         const HybridOCPDiscretization& discretization) {
  discretization.disp(os);
  return os;
}

} // namespace robotoc 