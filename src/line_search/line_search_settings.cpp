#include "robotoc/line_search/line_search_settings.hpp"

#include <cassert>
#include <stdexcept>
#include <iostream>


namespace robotoc {

LineSearchSettings::LineSearchSettings(const std::string& line_search_method,
                                       const double step_size_reduction_rate, 
                                       const double min_step_size, 
                                       const double armijo_control_rate,
                                       const double margin_rate,
                                       const double eps)
  : line_search_method(line_search_method),
    step_size_reduction_rate(step_size_reduction_rate),
    min_step_size(min_step_size),
    armijo_control_rate(armijo_control_rate),
    margin_rate(margin_rate),
    eps(eps) {
  try {
    if (line_search_method != "filter" && line_search_method != "merit-backtracking") {
      throw std::out_of_range("invalid value: line_search_method must be either \"filter\" or \"merit-backtracking\"!");
    }
    if (step_size_reduction_rate <= 0) {
      throw std::out_of_range("invalid value: step_size_reduction_rate must be positive!");
    }
    if (min_step_size <= 0) {
      throw std::out_of_range("invalid value: min_step_size must be positive!");
    }
    if (armijo_control_rate <= 0) {
      throw std::out_of_range("invalid value: min_step_size must be positive!");
    }
    if (margin_rate <= 0) {
      throw std::out_of_range("invalid value: margin_rate must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}

LineSearchSettings::LineSearchSettings() 
  : line_search_method(),
    step_size_reduction_rate(0),
    min_step_size(0),
    armijo_control_rate(0),
    margin_rate(0) {
}


LineSearchSettings::~LineSearchSettings() {
}

LineSearchSettings LineSearchSettings::defaultSettings() {
  LineSearchSettings s("filter", 0.75, 0.05, 0.001, 0.05, 1.0e-08);
  return s;
}

} // namespace robotoc