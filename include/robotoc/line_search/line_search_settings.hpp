#ifndef ROBOTOC_LINE_SEARCH_SETTINGS_HPP_
#define ROBOTOC_LINE_SEARCH_SETTINGS_HPP_

#include <iostream>


namespace robotoc {

/// 
/// @enum LineSearchMethod
/// @brief Type of the line search method.
///
enum class LineSearchMethod {
  Filter,
  MeritBacktracking,
};

///
/// @class LineSearchSettings
/// @brief Settings for the line search. 
///
struct LineSearchSettings {
public:
  ///
  /// @brief If set to LineSearchMethod::Filter, filter method is used as 
  /// a search scheme. If set to LineSearchMethod::MeritBacktracking, 
  /// backtracking line search is used.
  ///
  LineSearchMethod line_search_method = LineSearchMethod::Filter;

  ///
  /// @brief Reduction rate of the step size.
  ///
  double step_size_reduction_rate = 0.75;

  ///
  /// @brief Minimum step size.
  ///
  double min_step_size = 0.05;

  ///
  /// @brief The reduction rate of the cost in the filter line search method. 
  ///
  double filter_cost_reduction_rate = 0.005;

  ///
  /// @brief The reduction rate of the constraint violation in the filter line 
  /// search method. 
  ///
  double filter_constraint_violation_reduction_rate = 0.005;

  ///
  /// @brief Control rate in Armijo condition. This value 
  /// sets the slope of the linear approximation in the Armijo condition.
  ///
  double armijo_control_rate = 0.001;

  ///
  /// @brief Margin rate to determine penalty parameter in
  /// the merit function.
  ///
  double margin_rate = 0.05;

  ///
  /// @brief A small positive value. Used to calculate directional derivative.
  ///
  double eps = 1.0e-08;

  ///
  /// @brief Displays the line search settings onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const LineSearchSettings& line_search_settings);

};

} // namespace robotoc

#endif // ROBOTOC_LINE_SEARCH_SETTINGS_HPP_