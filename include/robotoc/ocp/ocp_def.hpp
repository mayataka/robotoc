#ifndef ROBOTOC_OCP_DEF_HPP_
#define ROBOTOC_OCP_DEF_HPP_

#include <memory>

#include "Eigen/Core"

#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"

namespace robotoc {

///
/// @class OCPDef
/// @brief A data structure for an optimal control problem.
///
struct OCPDef {
  std::shared_ptr<CostFunction> cost;
  std::shared_ptr<Constraints> constraints;
};

} // namespace robotoc

#endif // ROBOTOC_OCP_DEF_HPP_