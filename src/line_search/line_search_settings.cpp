#include "robotoc/line_search/line_search_settings.hpp"

#include <cassert>
#include <stdexcept>
#include <iostream>


namespace robotoc {

void LineSearchSettings::disp(std::ostream& os) const {
  os << "Line search settings:" << "\n";
  os << "  line search method: ";
  if (line_search_method == LineSearchMethod::Filter) os << "filter" << "\n";
  else os << "merit-backtracking" << "\n";
  os << "  step size reduction rate: " << step_size_reduction_rate << "\n";
  os << "  min step size: " << min_step_size << "\n";
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