#include "ocp/line_search_filter.hpp"


namespace idocp {

LineSearchFilter::LineSearchFilter() 
  : filter_() {
}

bool LineSearchFilter::isAccepted(const double cost, 
                                  const double constraint_residual) {
  for (auto pair : filter_) {
    if (cost > pair.first && constraint_residual > pair.second) {
      return false;
    }
  }
  return true;
}


void LineSearchFilter::append(const double cost, 
                              const double constraint_residual) {
  std::vector<std::pair<double, double>>::iterator it = filter_.begin();
  while (it != filter_.end()) {
    if (cost <= it->first && constraint_residual <= it->second) {
      it = filter_.erase(it);
    }
    else {
      ++it;
    }
  }
}


void LineSearchFilter::clear() {
  filter_.clear();
}

} // namespace idocp