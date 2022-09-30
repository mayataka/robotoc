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
  ///
  /// @brief Performance index.
  ///
  PerformanceIndex performance_index;

  ///
  /// @brief Cost function data.
  ///
  CostFunctionData cost_data;

  ///
  /// @brief Constraints data.
  ///
  ConstraintsData constraints_data;

  ///
  /// @brief State equation data.
  ///
  StateEquationData state_equation_data;

  ///
  /// @brief Contact dynamics data.
  ///
  ContactDynamicsData contact_dynamics_data;

  ///
  /// @brief Switching constraint data.
  ///
  SwitchingConstraintData switching_constraint_data;

  ///
  /// @brief Returns the lp norm of the primal feasibility, i.e., the constraint 
  /// violation. Default norm is l1-norm. You can also specify l-infty norm by 
  /// passing Eigen::Infinity as the template parameter.
  /// @tparam p Index of norm. Default is 1 (l1-norm).
  /// @return The lp norm of the primal feasibility.
  ///
  template <int p=1>
  double primalFeasibility() const {
    double feasibility = 0;
    feasibility += constraints_data.template primalFeasibility<p>();
    feasibility += contact_dynamics_data.template primalFeasibility<p>();
    return feasibility;
  }

  ///
  /// @brief Returns the lp norm of the dual feasibility. Default norm is 
  /// l1-norm. You can also specify l-infty norm by passing Eigen::Infinity as 
  /// the template parameter.
  /// @tparam p Index of norm. Default is 1 (l1-norm).
  /// @return The lp norm of the dual feasibility.
  ///
  template <int p=1>
  double dualFeasibility() const {
    double feasibility = 0;
    feasibility += constraints_data.template dualFeasibility<p>();
    feasibility += contact_dynamics_data.template dualFeasibility<p>();
    return feasibility;
  }

  ///
  /// @brief Returns the squared norm of the KKT residual, that is, 
  /// the primal and dual residual. 
  /// @return The squared norm of the KKT residual.
  ///
  double KKTError() const {
    double error = 0;
    error += contact_dynamics_data.KKTError();
    error += constraints_data.KKTError();
    return error;
  }

};

} // namespace robotoc

#endif // ROBOTOC_OCP_DATA_HPP_