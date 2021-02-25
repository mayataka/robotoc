#include "idocp/hybrid/parnmpc_discretizer.hpp"

#include <iostream>

namespace idocp {

void ParNMPCDiscretizer::showInfo() const {
  std::cout << "----- The discretized optimal control problem (ParNMPC) -----" << std::endl;
  std::cout << "T = " << T_ << std::endl;
  std::cout << "N_ideal = " << N_ideal() << std::endl;
  std::cout << "N = " << N() << std::endl;
  std::cout << "N_impulse = " << N_impulse() << std::endl;
  std::cout << "N_lift = " << N_lift() << std::endl;
  std::cout << "N_all = " << N_all() << std::endl;
  for (int i=0; i<N(); ++i) {
    if (isTimeStageAfterImpulse(i)) {
      const int index = impulseIndexBeforeTimeStage(i);
      std::cout << "impulse = " << index 
                << ": t =  " << t_impulse(index) 
                << ": dt_aux = " << dt_aux(index) << std::endl;
    }
    else if (isTimeStageAfterLift(i)) {
      const int index = liftIndexBeforeTimeStage(i);
      std::cout << "lift = " << index 
                << ": t =  " << t_lift(index) 
                << ": dt_lift = " << dt_lift(index) << std::endl;
    }
    std::cout << "i = " << i << ": t = " << t(i) << ": dt = " << dt(i) << std::endl; 
  }
}

} // namespace idocp 