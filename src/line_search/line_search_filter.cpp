#include "robotoc/line_search/line_search_filter.hpp"

#include <cassert>
#include <stdexcept>
#include <iostream>


namespace robotoc {

LineSearchFilter::LineSearchFilter(const double cost_reduction_rate, 
                                   const double constraints_reduction_rate) 
  : filter_(),
    cost_reduction_rate_(cost_reduction_rate),
    constraints_reduction_rate_(constraints_reduction_rate) {
  if (cost_reduction_rate_ <= 0) {
    throw std::out_of_range("[LineSearchFilter] invalid argument: cost_reduction_rate must be positive!");
  }
  if (constraints_reduction_rate_ <= 0) {
    throw std::out_of_range("[LineSearchFilter] invalid argument: constraints_reduction_rate must be positive!");
  }
}


LineSearchFilter::~LineSearchFilter() {
}


bool LineSearchFilter::isAccepted(const double cost, 
                                  const double constraint_violation) {
  assert(constraint_violation >= 0);
  if (filter_.empty()) {
    return true;
  }
  auto it = filter_.lower_bound(std::make_pair(cost, constraint_violation));
  if (it == filter_.begin()) {
    return true;
  }
  std::pair<double, double> prev = *(--it);
  if (cost < prev.first || constraint_violation < prev.second) {
    return true;
  }
  return false;
}


void LineSearchFilter::augment(const double cost, 
                               const double constraint_violation) {
  assert(constraint_violation >= 0);
  const double new_cost = cost - cost_reduction_rate_ * constraint_violation;
  const double new_constraint_violation = (1 - constraints_reduction_rate_) * constraint_violation;
  auto result = filter_.insert(std::make_pair(new_cost, new_constraint_violation));
  auto it = ++(result.first);
  while (it != filter_.end()) {
    if (cost <= it->first && constraint_violation <= it->second) {
        it = filter_.erase(it);
    }
    else {
      break;
    }
  }
}


void LineSearchFilter::clear() {
  filter_.clear();
}


bool LineSearchFilter::isEmpty() const {
  if (filter_.empty()) {
    return true;
  }
  else {
    return false;
  }
}

} // namespace robotoc