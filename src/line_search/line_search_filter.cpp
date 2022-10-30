#include "robotoc/line_search/line_search_filter.hpp"

#include <cassert>
#include <stdexcept>
#include <iostream>


namespace robotoc {

LineSearchFilter::LineSearchFilter(const double cost_reduction_rate, 
                                   const double constraint_violation_reduction_rate) 
  : filter_(),
    cost_reduction_rate_(cost_reduction_rate),
    constraint_violation_reduction_rate_(constraint_violation_reduction_rate) {
  if (cost_reduction_rate <= 0) {
    throw std::out_of_range("[LineSearchFilter] invalid argument: cost_reduction_rate must be positive!");
  }
  if (constraint_violation_reduction_rate <= 0) {
    throw std::out_of_range("[LineSearchFilter] invalid argument: constraint_violation_reduction_rate must be positive!");
  }
}


bool LineSearchFilter::isAccepted(const double cost, 
                                  const double constraint_violation) const {
  assert(constraint_violation >= 0);
  if (filter_.empty()) {
    return true;
  }
  for (const auto& e : filter_) {
    if ((cost < e.first - cost_reduction_rate_ * e.second) 
        || (constraint_violation < (1.0 - constraint_violation_reduction_rate_) * e.second)) {
      return true;
    }
  }
  return false;
}


void LineSearchFilter::augment(const double cost, 
                               const double constraint_violation) {
  assert(constraint_violation >= 0);
  if (!isAccepted(cost, constraint_violation)) {
    return;
  }

  auto it = filter_.begin();
  while (it != filter_.end()) {
    if ((it->first <= cost) && (it->second <= constraint_violation)) {
      it = filter_.erase(it);
    }
    else {
      ++it;
    }
  }

  filter_.emplace_back(cost, constraint_violation);
}


void LineSearchFilter::clear() {
  filter_.clear();
}


bool LineSearchFilter::isEmpty() const {
  return filter_.empty();
}

} // namespace robotoc