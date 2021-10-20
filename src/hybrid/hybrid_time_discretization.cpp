#include "idocp/hybrid/hybrid_ocp_discretization.hpp"


namespace idocp {

void HybridOCPDiscretization::disp(std::ostream& os) const {
  os << "The discretized optimal control problem (OCP):" << std::endl;
  os << "  T = " << T_ << std::endl;
  os << "  N_ideal = " << N_ideal() << std::endl;
  os << "  N = " << N() << std::endl;
  os << "  N_impulse = " << N_impulse() << std::endl;
  os << "  N_lift = " << N_lift() << std::endl;
  os << "  N_all = " << N_all() << std::endl;
  os << "  isFormulationTractable: ";
  if (isFormulationTractable()) os << "true" << std::endl;
  else os << "false" << std::endl;
  os << "  isSwitchingTimeConsistent: ";
  if (isSwitchingTimeConsistent()) os << "true" << std::endl;
  else os << "false" << std::endl;
  for (int i=0; i<N(); ++i) {
    os << "  i = " << i << ": t = " << t(i) << ": dt = " << dt(i) << std::endl; 
    if (isTimeStageBeforeImpulse(i)) {
      const int index = impulseIndexAfterTimeStage(i);
      os << "  impulse = " << index 
                << ": t =  " << t_impulse(index) 
                << ": dt_aux = " << dt_aux(index) << std::endl;
    }
    else if (isTimeStageBeforeLift(i)) {
      const int index = liftIndexAfterTimeStage(i);
      os << "  lift = " << index 
                << ": t =  " << t_lift(index) 
                << ": dt_lift = " << dt_lift(index) << std::endl;
    }
  }
  os << "  i = " << N() << ": t = " << t(N()) << std::flush; 
}


std::ostream& operator<<(std::ostream& os, 
                         const HybridOCPDiscretization& discretization) {
  discretization.disp(os);
  return os;
}

} // namespace idocp 