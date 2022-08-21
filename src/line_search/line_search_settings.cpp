#include "robotoc/line_search/line_search_settings.hpp"

#include <cassert>
#include <stdexcept>
#include <iostream>


namespace robotoc {

LineSearchSettings::LineSearchSettings(const LineSearchMethod _line_search_method,
                                       const double _step_size_reduction_rate, 
                                       const double _min_step_size, 
                                       const double _armijo_control_rate,
                                       const double _margin_rate,
                                       const double _eps)
  : line_search_method(_line_search_method),
    step_size_reduction_rate(_step_size_reduction_rate),
    min_step_size(_min_step_size),
    armijo_control_rate(_armijo_control_rate),
    margin_rate(_margin_rate),
    eps(_eps) {
  if (step_size_reduction_rate <= 0) {
    throw std::out_of_range("[LineSearchSettings] invalid argument: 'step_size_reduction_rate' must be positive!");
  }
  if (min_step_size <= 0) {
    throw std::out_of_range("[LineSearchSettings] invalid argument: 'min_step_size' must be positive!");
  }
  if (armijo_control_rate <= 0) {
    throw std::out_of_range("[LineSearchSettings] invalid argument: 'min_step_size' must be positive!");
  }
  if (margin_rate <= 0) {
    throw std::out_of_range("[LineSearchSettings] invalid argument: 'margin_rate' must be positive!");
  }
}



LineSearchSettings::LineSearchSettings() 
  : LineSearchSettings(LineSearchMethod::Filter, 0.75, 0.05, 0.001, 0.05, 1.0e-08) {
}


LineSearchSettings::~LineSearchSettings() {
}


LineSearchSettings LineSearchSettings::defaultSettings() {
  return LineSearchSettings();
}


void LineSearchSettings::disp(std::ostream& os) const {
  os << "Line search settings:" << std::endl;
  os << "  line search method: ";
  if (line_search_method == LineSearchMethod::Filter) os << "filter" << std::endl;
  else os << "merit-backtracking" << std::endl;
  os << "  step size reduction rate: " << step_size_reduction_rate << std::endl;
  os << "  min step size: " << min_step_size << std::endl;
  os << "  armijo control rate: " << armijo_control_rate << std::endl;
  os << "  margin rate: " << margin_rate << std::endl;
  os << "  eps: " << eps << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const LineSearchSettings& line_search_settings) {
  line_search_settings.disp(os);
  return os;
}

} // namespace robotoc