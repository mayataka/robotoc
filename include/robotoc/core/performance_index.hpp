#ifndef ROBOTOC_PERFORMANCE_INDEX_HPP_
#define ROBOTOC_PERFORMANCE_INDEX_HPP_

#include <iostream>


namespace robotoc {

///
/// @typedef PerformanceIndex
/// @brief Performance index of optimal control problems.
///
struct PerformanceIndex {
  ///
  /// @brief Cost value.
  ///
  double cost = 0.0;

  ///
  /// @brief Cost value of slack barrier function.
  ///
  double cost_barrier = 0.0;

  ///
  /// @brief Primal feasiblity.
  ///
  double primal_feasibility = 0.0;

  ///
  /// @brief Dual feasiblity.
  ///
  double dual_feasibility = 0.0;

  ///
  /// @brief KKT error.
  ///
  double kkt_error = 0.0;

  ///
  /// @brief Overload the operator.
  ///
  PerformanceIndex& operator+=(const PerformanceIndex& other);

  ///
  /// @brief Overload the operator.
  ///
  PerformanceIndex operator+(const PerformanceIndex& other) const;

  ///
  /// @brief Displays the performance index onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const PerformanceIndex& kkt_residual);
};

} // namespace robotoc

#endif // ROBOTOC_PERFORMANCE_INDEX_HPP_