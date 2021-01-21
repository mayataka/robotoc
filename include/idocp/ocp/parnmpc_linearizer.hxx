#ifndef IDOCP_PARNMPC_LINEARIZER_HXX_ 
#define IDOCP_PARNMPC_LINEARIZER_HXX_

#include "idocp/ocp/parnmpc_linearizer.hpp"

#include <cassert>


namespace idocp {

inline const Eigen::VectorXd& ParNMPCLinearizer::q_prev(
    const ParNMPCDiscretizer& discretizer, const Eigen::VectorXd& q, 
    const Solution& s, const int time_stage) {
  assert(time_stage >= 0);
  assert(time_stage < discretizer.N());
  if (discretizer.isTimeStageAfterImpulse(time_stage)) {
    return s.impulse[discretizer.impulseIndexBeforeTimeStage(time_stage)].q;
  }
  else if (discretizer.isTimeStageAfterLift(time_stage)) {
    return s.lift[discretizer.liftIndexBeforeTimeStage(time_stage)].q;
  }
  else if (time_stage > 0) {
    return s[time_stage-1].q;
  }
  else { 
    assert(time_stage == 0);
    return q;
  }
}


inline const Eigen::VectorXd& ParNMPCLinearizer::v_prev(
    const ParNMPCDiscretizer& discretizer, const Eigen::VectorXd& v, 
    const Solution& s, const int time_stage) {
  assert(time_stage >= 0);
  assert(time_stage < discretizer.N());
  if (discretizer.isTimeStageAfterImpulse(time_stage)) {
    return s.impulse[discretizer.impulseIndexBeforeTimeStage(time_stage)].v;
  }
  else if (discretizer.isTimeStageAfterLift(time_stage)) {
    return s.lift[discretizer.liftIndexBeforeTimeStage(time_stage)].v;
  }
  else if (time_stage > 0) {
    return s[time_stage-1].v;
  }
  else { 
    assert(time_stage == 0);
    return v;
  }
}

} // namespace idocp 

#endif // IDOCP_PARNMPC_LINEARIZER_HXX_ 