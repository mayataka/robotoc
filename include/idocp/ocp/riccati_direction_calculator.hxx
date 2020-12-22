#ifndef IDOCP_RICCATI_DIRECTION_CALCULATOR_HXX_
#define IDOCP_RICCATI_DIRECTION_CALCULATOR_HXX_

#include "idocp/ocp/riccati_direction_calculator.hpp"

namespace idocp {

template <>
inline void RiccatiDirectionCalculator::computeInitialStateDirection<true>(
    const std::vector<Robot>& robots, const Eigen::VectorXd& q, 
    const Eigen::VectorXd& v, const Solution& s, Direction& d, 
    const double sampling_period) {
  assert(q.size() == robots[0].dimq());
  assert(v.size() == robots[0].dimv());
  assert(sampling_period > 0);
  robots[0].subtractConfiguration(q, s[0].q, d[0].dq());
  d[0].dv() = v - s[0].v;
  d[0].dq().noalias() += sampling_period * v;
  d[0].dv().noalias() += sampling_period * s[0].a;
}


template <>
inline void RiccatiDirectionCalculator::computeInitialStateDirection<false>(
    const std::vector<Robot>& robots, const Eigen::VectorXd& q, 
    const Eigen::VectorXd& v, const Solution& s, Direction& d, 
    const double sampling_period) {
  assert(q.size() == robots[0].dimq());
  assert(v.size() == robots[0].dimv());
  robots[0].subtractConfiguration(q, s[0].q, d[0].dq());
  d[0].dv() = v - s[0].v;
}


inline const SplitRiccatiFactorization& RiccatiDirectionCalculator::
next_riccati_factorization(const OCPDiscretizer& ocp_discretizer, 
                           const RiccatiFactorization& factorization, 
                           const int time_stage) {
  assert(time_stage >= 0);
  assert(time_stage < ocp_discretizer.N());
  if (ocp_discretizer.isTimeStageBeforeImpulse(time_stage)) {
    return factorization.impulse[ocp_discretizer.impulseIndexAfterTimeStage(time_stage)];
  }
  else if (ocp_discretizer.isTimeStageBeforeLift(time_stage)) {
    return factorization.lift[ocp_discretizer.liftIndexAfterTimeStage(time_stage)];
  }
  else {
    return factorization[time_stage+1];
  }
}

} // namespace idocp 

#endif // IDOCP_RICCATI_DIRECTION_CALCULATOR_HXX_ 