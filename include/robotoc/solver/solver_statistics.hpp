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
class SolverStatistics {
public:
  ///
  /// @brief Default constructor. 
  ///
  SolverStatistics();

  ///
  /// @brief Destructor. 
  ///
  ~SolverStatistics();

  ///
  /// @brief Default copy constructor. 
  ///
  SolverStatistics(const SolverStatistics&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  SolverStatistics& operator=(const SolverStatistics&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SolverStatistics(SolverStatistics&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SolverStatistics& operator=(SolverStatistics&&) noexcept = default;

  ///
  /// @brief Flags whether the convergence is achieved.
  ///
  bool convergence;

  ///
  /// @brief Number of iterations until convergence.
  ///
  int iter;

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