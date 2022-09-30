#ifndef ROBOTOC_OCP_DATA_HPP_
#define ROBOTOC_OCP_DATA_HPP_

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
/// @class OCPData
/// @brief A data structure for an optimal control problem.
///
struct OCPData {

  PerformanceIndex performance_index;

  CostFunctionData cost_data;

  ConstraintsData constraints_data;

  StateEquationData state_equation_data;

  ContactDynamicsData contact_dynamics_data;

  SwitchingConstraintData switching_constraint_data;

  template <int p=1>
  double primalFeasibility() const {
    double feasibility = 0;
    feasibility += constraints_data.template primalFeasibility<p>();
    feasibility += contact_dynamics_data.template primalFeasibility<p>();
    return feasibility;
  }

  template <int p=1>
  double dualFeasibility() const {
    double feasibility = 0;
    feasibility += constraints_data.template dualFeasibility<p>();
    feasibility += contact_dynamics_data.template dualFeasibility<p>();
    return feasibility;
  }

  double KKTError() const {
    double error = 0;
    std::cout << "contact dynamics KKT: " << contact_dynamics_data.KKTError() << std::endl;
    std::cout << "constraints KKT: " << constraints_data.KKTError() << std::endl;
    error += contact_dynamics_data.KKTError();
    error += constraints_data.KKTError();
    return error;
  }

};

} // namespace robotoc

#endif // ROBOTOC_OCP_DATA_HPP_