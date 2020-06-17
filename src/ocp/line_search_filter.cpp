#include "ocp/line_search_filter.hpp"


namespace idocp {

LineSearchFilter::LineSearchFilter() 
  : filter_(),
    cost_reduction_rate_(0.005),
    constraints_reduction_rate_(0.005) {
}

bool LineSearchFilter::isAccepted(const double cost, 
                                  const double constraint_violation) {
  if (!filter_.empty()) {
    for (auto pair : filter_) {
      if (cost > pair.first && constraint_violation > pair.second) {
        return false;
      }
    }
  }
  return true;
}


void LineSearchFilter::append(const double cost, 
                              const double constraint_violation) {
  if (!filter_.empty()) {
    std::vector<std::pair<double, double>>::iterator it = filter_.begin();
    while (it != filter_.end()) {
      if (cost <= it->first && constraint_violation <= it->second) {
        it = filter_.erase(it);
      }
      else {
        ++it;
      }
    }
  }
  filter_.push_back(std::make_pair(
      cost-cost_reduction_rate_*constraint_violation, 
      (1-constraints_reduction_rate_)*constraint_violation));
}


void LineSearchFilter::clear() {
  filter_.clear();
}

} // namespace idocp