#include "idocp/hybrid/hybrid_time_discretization.hpp"

#include <iostream>

namespace idocp {

void HybridTimeDiscretization::showInfo() const {
  std::cout << "----- The discretized optimal control problem (OCP) -----" << std::endl;
  std::cout << "T = " << T_ << std::endl;
  std::cout << "N_ideal = " << N_ideal() << std::endl;
  std::cout << "N = " << N() << std::endl;
  std::cout << "N_impulse = " << N_impulse() << std::endl;
  std::cout << "N_lift = " << N_lift() << std::endl;
  std::cout << "N_all = " << N_all() << std::endl;
  std::cout << "isFormulationTractable: ";
  if (isFormulationTractable()) std::cout << "true" << std::endl;
  else std::cout << "false" << std::endl;
  std::cout << "isSwitchingTimeConsistent: ";
  if (isSwitchingTimeConsistent()) std::cout << "true" << std::endl;
  else std::cout << "false" << std::endl;
  for (int i=0; i<N(); ++i) {
    std::cout << "i = " << i << ": t = " << t(i) << ": dt = " << dt(i) << std::endl; 
    if (isTimeStageBeforeImpulse(i)) {
      const int index = impulseIndexAfterTimeStage(i);
      std::cout << "impulse = " << index 
                << ": t =  " << t_impulse(index) 
                << ": dt_aux = " << dt_aux(index) << std::endl;
    }
    else if (isTimeStageBeforeLift(i)) {
      const int index = liftIndexAfterTimeStage(i);
      std::cout << "lift = " << index 
                << ": t =  " << t_lift(index) 
                << ": dt_lift = " << dt_lift(index) << std::endl;
    }
  }
  std::cout << "i = " << N() << ": t = " << t(N()) << std::endl; 
  std::cout << "---------------------------------------------------------" << std::endl;
}

} // namespace idocp 