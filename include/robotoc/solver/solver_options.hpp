#ifndef ROBOTOC_SOLVER_OPTIONS_HPP_
#define ROBOTOC_SOLVER_OPTIONS_HPP_

#include <iostream>

#include "robotoc/ocp/discretization_method.hpp"
#include "robotoc/line_search/line_search_settings.hpp"
#include "robotoc/solver/interpolation_order.hpp"


namespace robotoc {

///
/// @class SolverOptions
/// @brief Options of optimal control solvers. 
///
struct SolverOptions {
  ///
  /// @brief Maximum number of iterations. Must be non-negative. 
  /// Default is 100. 
  ///
  int max_iter = 100;

  ///
  /// @brief Tolerance of the l2-norm of the (perturbed) KKT residual for 
  /// convergence. Default is 1.0e-07. Must be positive.
  ///
  double kkt_tol = 1.0e-07;

  ///
  /// @brief Initial barrier parameter. Must be positive. Default is 1.0e-03.
  ///
  double mu_init = 1.0e-03;

  ///
  /// @brief Minimum barrier parameter. Must be positive. Default is 1.0e-03
  /// (that is, no barrier parameter reduction with the default settings).
  ///
  double mu_min = 1.0e-03;

  ///
  /// @brief Tolerance of the l2-norm of the (perturbed) KKT residual for 
  /// the barrier parameter update. Default is 1.0e-07. Must be positive.
  ///
  double kkt_tol_mu = 1.0e-07;

  ///
  /// @brief Linear decrease rate of barrier parameter. Default is 0.2.
  /// @note New barrier parameter mu is obtained by taking the minimum of 
  /// mu*"mu_linear_decrease_factor" and mu^"superlinear_decrease_power".
  ///
  double mu_linear_decrease_factor = 0.2;

  ///
  /// @brief Superlinear decrease rate of barrier parameter. Default is 1.5.
  /// @note New barrier parameter mu is obtained by taking the minimum of 
  /// mu*"mu_linear_decrease_factor" and mu^"superlinear_decrease_power".
  ///
  double mu_superlinear_decrease_power = 1.5;

  ///
  /// @brief Flag to enable the line search. Default is false.
  ///
  bool enable_line_search = false;

  ///
  /// @brief Line search settings. 
  ///
  LineSearchSettings line_search_settings;

  ///
  /// @brief Discretization method of the hybrid optimal control problem.
  /// Only used in OCPSolver without the STO problem. Default is 
  /// DiscretizationMethod::GridBased.
  /// @note For the STO problem, discretization method is fixed to 
  /// DiscretizationMethod::PhaseBased regardless of this value.
  ///
  DiscretizationMethod discretization_method = DiscretizationMethod::GridBased;

  ///
  /// @brief Number of initial inner iterations in which a large regularization 
  /// for the STO problem is added, where the inner iteration means the 
  /// iterations after each mesh-refinement. Default is 0
  /// (that is, there is no iterations with large STO regularization).
  /// @note The purpose of a large regularization on the STO problem is to  
  /// stabilize the numerical computation. With sufficiently large  
  /// regularization on the STO, the optimization problem involving the STO 
  /// problem is the same as one without STO problem.
  ///
  int initial_sto_reg_iter = 0;

  ///
  /// @brief Magnitude of the large regularization for the STO problem at the 
  /// initial iterations speficied by initial_sto_reg_iter. 
  /// Default is 1.0e30. 
  /// @note The purpose of a large regularization on the STO problem is to  
  /// stabilize the numerical computation. With sufficiently large  
  /// regularization on the STO, the optimization problem involving the STO 
  /// problem is the same as one without STO problem.
  ///
  double initial_sto_reg = 1.0e30;

  ///
  /// @brief Tolerance of the KKT residual to perform mesh-refinement in the
  /// STO problem. Default is 1.0e-01. Must be positive.
  ///
  double kkt_tol_mesh = 0.1;

  ///
  /// @brief Tolerance of the maximum discretization step size to perform  
  /// mesh-refinement in the STO problem. Default is 0.
  ///
  double max_dt_mesh = 0.0;

  ///
  /// @brief Maximum magnitude of the nominal direction of the switching time 
  /// computed by the Riccati recursion algorithm. Used in a heuristic 
  /// regularization on the dynamic programming recursion, which is similar to
  /// the trust-region strategy (this value is a kind of a trust region w.r.t. 
  /// the switching time direction). Must be positive. 
  /// The regularization increases as this value is closer to 0. A larger 
  /// regularization means that the numerical computation is more stable  
  /// while the convergence speed can become slower. Default is 0.1, that is,
  /// a regularization is added so that the magnitude of the nominal switching  
  /// time direction is limited to under 0.1 in default.
  ///
  double max_dts_riccati = 0.1;

  ///
  /// @brief If true, the solution initial guess is constructed from the linear 
  /// interpolation of the previous solution. 
  ///
  bool enable_solution_interpolation = true;

  ///
  /// @brief Order of the solution interpolation if enable_solution_interpolation
  /// is true. 
  ///
  InterpolationOrder interpolation_order = InterpolationOrder::Linear;

  ///
  /// @brief If true, the CPU time is measured at each solve().
  ///
  bool enable_benchmark = false;

  ///
  /// @brief Displays the solver settings onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const SolverOptions& solver_options);

};

} // namespace robotoc

#endif // ROBOTOC_SOLVER_OPTIONS_HPP_