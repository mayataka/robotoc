#ifndef IDOCP_PARNMPC_LINEARIZER_HXX_ 
#define IDOCP_PARNMPC_LINEARIZER_HXX_

#include "idocp/ocp/parnmpc_linearizer.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>


namespace idocp {

inline const Eigen::VectorXd& ParNMPCLinearizer::q_prev(
    const ParNMPCDiscretizer& discretizer, const Eigen::VectorXd& q, 
    const Solution& s, const int time_stage) {
  assert(time_stage >= 0);
  assert(time_stage <= discretizer.N());
  if (time_stage == 0) {
    return q;
  }
  else if (discretizer.isTimeStageBeforeImpulse(time_stage-1)) {
    return s.aux[discretizer.impulseIndexAfterTimeStage(time_stage-1)].q;
  }
  else if (discretizer.isTimeStageBeforeLift(time_stage-1)) {
    return s.lift[discretizer.liftIndexAfterTimeStage(time_stage-1)].q;
  }
  else {
    return s[time_stage-1].q;
  }
}


inline const Eigen::VectorXd& ParNMPCLinearizer::v_prev(
    const ParNMPCDiscretizer& discretizer, const Eigen::VectorXd& v, 
    const Solution& s, const int time_stage) {
  assert(time_stage >= 0);
  assert(time_stage <= discretizer.N());
  if (time_stage == 0) {
    return v;
  }
  else if (discretizer.isTimeStageBeforeImpulse(time_stage-1)) {
    return s.aux[discretizer.impulseIndexAfterTimeStage(time_stage-1)].v;
  }
  else if (discretizer.isTimeStageBeforeLift(time_stage-1)) {
    return s.lift[discretizer.liftIndexAfterTimeStage(time_stage-1)].v;
  }
  else {
    return s[time_stage-1].v;
  }
}

} // namespace idocp 

#endif // IDOCP_PARNMPC_LINEARIZER_HXX_ 