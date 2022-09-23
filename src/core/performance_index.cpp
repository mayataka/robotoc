#include "robotoc/core/performance_index.hpp"

#include <iostream>


namespace robotoc {

void PerformanceIndex::disp(std::ostream& os) const {
  os << "PerformanceIndex:" << std::endl;
  os << "  cost =             = " << cost << std::endl;
  os << "  cost_barrier       = " << cost_barrier << std::endl;
  os << "  primal_feasibility = " << primal_feasibility << std::endl;
  os << "  dual_feasibility   = " << dual_feasibility << std::endl;
  os << "  kkt_error          = " << kkt_error << std::flush;
}


std::ostream& operator<<(std::ostream& os, const PerformanceIndex& peformance) {
  peformance.disp(os);
  return os;
}


} // namespace robotoc