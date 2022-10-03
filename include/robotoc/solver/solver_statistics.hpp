#ifndef ROBOTOC_SOLVER_STATISTICS_HPP_
#define ROBOTOC_SOLVER_STATISTICS_HPP_

#include <vector>
#include <deque>
#include <iostream>


namespace robotoc {

///
/// @class SolverStatistics
/// @brief Statistics of optimal control solvers. 
///
struct SolverStatistics {
  ///
  /// @brief Flags whether the convergence is achieved.
  ///
  bool convergence = false;

  ///
  /// @brief Number of iterations until convergence.
  ///
  int iter = 0;

  ///
  /// @brief l2-norm of the KKT residual at each iteration.
  ///
  std::vector<double> kkt_error;

  ///
  /// @brief Primal step sizes at each iteration.
  ///
  std::vector<double> primal_step_size;

  ///
  /// @brief Dual step sizes at each iteration.
  ///
  std::vector<double> dual_step_size;

  ///
  /// @brief Switching times at each iteration.
  ///
  std::vector<std::deque<double>> ts;

  ///
  /// @brief Iterations where the mesh-refinements are carried out.
  ///
  std::vector<int> mesh_refinement_iter;

  ///
  /// @brief CPU time is stored if SolverOptions::enable_benchmark is true.
  ///
  double cpu_time = 0;

  ///
  /// @brief Reserves the data.
  /// @param[in] size Size of the new data.
  ///
  void reserve(const int size);

  ///
  /// @brief Clear the all elements.
  ///
  void clear();

  ///
  /// @brief Displays the solver statistics onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const SolverStatistics& solver_statistics);

};

} // namespace robotoc

#endif // ROBOTOC_SOLVER_STATISTICS_HPP_ 