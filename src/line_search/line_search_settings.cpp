#include "robotoc/line_search/line_search_settings.hpp"

#include <cassert>
#include <stdexcept>
#include <iostream>


namespace robotoc {

void LineSearchSettings::disp(std::ostream& os) const {
  auto lineSearchMethodToString = [](const LineSearchMethod& method) {
    switch (method)
    {
    case LineSearchMethod::Filter:
      return "Filter";
      break;
    case LineSearchMethod::MeritBacktracking:
      return "MeritBacktracking";
      break;
    default:
      return "";
      break;
    }
  };
  os << "Line search settings:" << "\n";
  os << "  line search method: " << lineSearchMethodToString(line_search_method) << "\n";
  os << "  step size reduction rate: " << step_size_reduction_rate << "\n";
  os << "  min step size: " << min_step_size << "\n";
  os << "  filter cost reduction rate: " << filter_cost_reduction_rate << "\n";
  os << "  filter constraint violation reduction rate: " << filter_constraint_violation_reduction_rate << "\n";
  os << "  armijo control rate: " << armijo_control_rate << "\n";
  os << "  margin rate: " << margin_rate << "\n";
  os << "  eps: " << eps << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const LineSearchSettings& line_search_settings) {
  line_search_settings.disp(os);
  return os;
}

} // namespace robotoc