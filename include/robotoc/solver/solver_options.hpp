#ifndef ROBOTOC_SOLVER_OPTIONS_HPP_
#define ROBOTOC_SOLVER_OPTIONS_HPP_

#include <iostream>

#include "robotoc/hybrid/discretization_method.hpp"
#include "robotoc/line_search/line_search_settings.hpp"


namespace robotoc {

///
/// @class SolverOptions
/// @brief Options of optimal control solvers. 
///
class SolverOptions {
public:
  ///
  /// @brief Constructor. 
  ///
  SolverOptions(const int max_iter, const double kkt_tol, const double mu_init,
                const double mu_min, const double kkt_tol_mu, 
                const double mu_linear_decrease_factor, 
                const double mu_superlinear_decrease_power, 
                const bool enable_line_search, 
                const LineSearchSettings& line_search_settings, 
                const DiscretizationMethod& discretization_method, 
                const int initial_sto_reg_iter, const double initial_sto_reg, 
                const double kkt_tol_mesh, const double max_dt_mesh, 
                const double max_dts_riccati, const int print_level);

  ///
  /// @brief Default constructor. 
  ///
  SolverOptions();

  ///
  /// @brief Destructor. 
  ///
  ~SolverOptions();

  ///
  /// @brief Default copy constructor. 
  ///
  SolverOptions(const SolverOptions&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  SolverOptions& operator=(const SolverOptions&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SolverOptions(SolverOptions&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SolverOptions& operator=(SolverOptions&&) noexcept = default;

  ///
  /// @brief Maximum number of iterations. Must be non-negative. 
  /// Default is 100. 
  ///
  int max_iter;

  ///
  /// @brief Tolerance of the l2-norm of the (perturbed) KKT residual for 
  /// convergence. Default is 1.0e-07. Must be positive.
  ///
  double kkt_tol;

  ///
  /// @brief Initial barrier parameter. Must be positive. Default is 1.0e-03.
  ///
  double mu_init;

  ///
  /// @brief Minimum barrier parameter. Must be positive. Default is 1.0e-03
  /// (that is, no barrier parameter reduction with the default settings).
  ///
  double mu_min;

  ///
  /// @brief Tolerance of the l2-norm of the (perturbed) KKT residual for 
  /// the barrier parameter update. Default is 1.0e-07. Must be positive.
  ///
  double kkt_tol_mu;

  ///
  /// @brief Linear decrease rate of barrier parameter. Default is 0.2.
  /// @note New barrier parameter mu is obtained by taking the minimum of 
  /// mu*"mu_linear_decrease_factor" and mu^"superlinear_decrease_power".
  ///
  double mu_linear_decrease_factor;

  ///
  /// @brief Superlinear decrease rate of barrier parameter. Default is 1.5.
  /// @note New barrier parameter mu is obtained by taking the minimum of 
  /// mu*"mu_linear_decrease_factor" and mu^"superlinear_decrease_power".
  ///
  double mu_superlinear_decrease_power;

  ///
  /// @brief Flag to enable the line search. Default is false.
  ///
  bool enable_line_search;

  ///
  /// @brief Line search settings. Default is 
  /// LineSearchSettings::defaultSettings().
  ///
  LineSearchSettings line_search_settings;

  ///
  /// @brief Discretization method of the hybrid optimal control problem.
  /// Only used in OCPSolver without the STO problem. Default is 
  /// DiscretizationMethod::GridBased.
  /// @note For the STO problem, discretization method is fixed to 
  /// DiscretizationMethod::PhaseBased regardless of this value.
  ///
  DiscretizationMethod discretization_method;

  ///
  /// @brief Number of initial iterations in which a large regularization for 
  /// the STO problem is added. Default is 0. (there is no iterations 
  /// with large STO regularization).
  /// @note The purpose of a large regularization on the STO problem is to  
  /// stabilize the numerical computation. With sufficiently large  
  /// regularization on the STO, the optimization problem involving the STO 
  /// problem is the same as one without STO problem.
  ///
  int initial_sto_reg_iter;

  ///
  /// @brief Magnitude of the large regularization for the STO problem at the 
  /// initial iterations speficied by initial_sto_reg_iter. 
  /// Default is 1.0e30. 
  /// @note The purpose of a large regularization on the STO problem is to  
  /// stabilize the numerical computation. With sufficiently large  
  /// regularization on the STO, the optimization problem involving the STO 
  /// problem is the same as one without STO problem.
  ///
  double initial_sto_reg;

  ///
  /// @brief Tolerance of the KKT residual to perform mesh-refinement in the
  /// STO problem. Default is 1.0e-01. Must be positive.
  ///
  double kkt_tol_mesh;

  ///
  /// @brief Tolerance of the maximum discretization step size to perform  
  /// mesh-refinement in the STO problem. Default is 0.
  ///
  double max_dt_mesh;

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
  double max_dts_riccati;

  ///
  /// @brief Print level. 
  /// 0: do not print anything. 
  /// 1: print the KKT residual at each iterations.
  /// 2: further prints the step sizes, etc.
  /// Default is 1.
  ///
  int print_level;

  ///
  /// @brief Returns options with default parameters.
  ///
  static SolverOptions defaultOptions();

  ///
  /// @brief Displays the solver settings onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const SolverOptions& solver_options);

};

} // namespace robotoc

#endif // ROBOTOC_SOLVER_OPTIONS_HPP_