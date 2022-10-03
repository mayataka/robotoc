#ifndef ROBOTOC_UNCONSTR_OCP_DATA_HPP_
#define ROBOTOC_UNCONSTR_OCP_DATA_HPP_

#include "robotoc/core/performance_index.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/constraints/constraints_data.hpp"
#include "robotoc/dynamics/unconstr_dynamics.hpp"


namespace robotoc {

///
/// @class UnconstrOCPData
/// @brief Data structure for the optimal control problem of unconstrained 
/// rigid-body systems.
///
struct UnconstrOCPData {
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
  UnconstrDynamics unconstr_dynamics;

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
    feasibility += unconstr_dynamics.template primalFeasibility<p>();
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
    return feasibility;
  }

  ///
  /// @brief Returns the squared norm of the KKT residual, that is, 
  /// the primal and dual residual. 
  /// @return The squared norm of the KKT residual.
  ///
  double KKTError() const {
    double error = 0;
    error += constraints_data.KKTError();
    error += unconstr_dynamics.KKTError();
    return error;
  }
};

} // namespace robotoc

#endif // ROBOTOC_UNCONSTR_OCP_DATA_HPP_