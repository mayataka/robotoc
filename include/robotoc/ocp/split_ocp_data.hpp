#ifndef ROBOTOC_SPLIT_OCP_DATA_HPP_
#define ROBOTOC_SPLIT_OCP_DATA_HPP_

#include <memory>

#include "Eigen/Core"

#include "robotoc/core/performance_index.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/constraints/constraints_data.hpp"
#include "robotoc/dynamics/state_equation_data.hpp"
#include "robotoc/dynamics/contact_dynamics_data.hpp"
#include "robotoc/dynamics/switching_constraint_data.hpp"


namespace robotoc {

///
/// @class SplitOCPData
/// @brief A data structure for an optimal control problem.
///
struct SplitOCPData {
  PerformanceIndex performance_index;
  CostFunctionData cost_data;
  ConstraintsData constraints_data;
  StateEquationData state_equation_data;
  ContactDynamicsData contact_dynamics_data;
  SwitchingConstraintData switching_constraint_data;
};

} // namespace robotoc

#endif // ROBOTOC_SPLIT_OCP_DATA_HPP_