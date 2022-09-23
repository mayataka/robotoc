#include "robotoc/core/performance_index.hpp"

#include <iostream>


namespace robotoc {

PerformanceIndex& PerformanceIndex::operator+=(const PerformanceIndex& other) {
  this->cost               += other.cost;
  this->cost_barrier       += other.cost_barrier;
  this->primal_feasibility += other.primal_feasibility;
  this->dual_feasibility   += other.dual_feasibility;
  this->kkt_error          += other.kkt_error;
  return *this;
}


PerformanceIndex PerformanceIndex::operator+(const PerformanceIndex& other) const {
  PerformanceIndex ret = *this;
  ret += other;
  return ret;
}


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