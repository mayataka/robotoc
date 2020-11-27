#ifndef IDOCP_OCP_SOLUTION_INTEGRATOR_HXX
#define IDOCP_OCP_SOLUTION_INTEGRATOR_HXX_

#include "idocp/ocp/ocp_solution_integrator.hpp"

namespace idocp {

inline bool OCPSolutionIntegrator::is_state_constraint_valid(
    const int time_stage_before_impulse) {
  assert(time_stage_before_impulse >= 0);
  if (time_stage_before_impulse > 0) {
    return true;
  }
  else {
    return false;
  }
}

} // namespace idocp

#endif // IDOCP_OCP_SOLUTION_INTEGRATOR_HXX_ 