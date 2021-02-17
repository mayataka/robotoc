#ifndef IDOCP_STATE_CONSTRAINT_JACOBIAN_HPP_ 
#define IDOCP_STATE_CONSTRAINT_JACOBIAN_HPP_

#include <vector>

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_state_constraint_jacobian.hpp"

namespace idocp {

class StateConstraintJacobian {
public: 
  StateConstraintJacobian(const Robot& robot, const int N_impulse=0) 
    : data(N_impulse, SplitStateConstraintJacobian(robot)) {
  }

  StateConstraintJacobian() 
    : data() {
  }

  StateConstraintJacobian(const StateConstraintJacobian&) = default;

  StateConstraintJacobian& operator=(const StateConstraintJacobian&) = default;

  StateConstraintJacobian(StateConstraintJacobian&&) noexcept = default;

  StateConstraintJacobian& operator=(StateConstraintJacobian&&) noexcept = default;

  SplitStateConstraintJacobian& operator[] (const int i) {
    assert(i >= 0);
    assert(i < data.size());
    return data[i];
  }

  const SplitStateConstraintJacobian& operator[] (const int i) const {
    assert(i >= 0);
    assert(i < data.size());
    return data[i];
  }

  std::vector<SplitStateConstraintJacobian> data;
};

} // namespace idocp

#endif // IDOCP_STATE_CONSTRAINT_JACOBIAN_HPP_