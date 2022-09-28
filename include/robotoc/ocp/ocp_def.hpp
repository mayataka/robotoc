#ifndef ROBOTOC_OCP_DEF_HPP_
#define ROBOTOC_OCP_DEF_HPP_

#include <memory>

#include "Eigen/Core"

#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/planner/contact_sequence.hpp"

namespace robotoc {

///
/// @class OCPDef
/// @brief A definition of an optimal control problem.
///
struct OCPDef {
  Robot robot;
  std::shared_ptr<CostFunction> cost;
  std::shared_ptr<Constraints> constraints;
  std::shared_ptr<ContactSequence> contact_sequence;
  double T;
  int N;
};

} // namespace robotoc

#endif // ROBOTOC_OCP_DEF_HPP_