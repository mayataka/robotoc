#include "robotoc/line_search/line_search_filter.hpp"

#include <cassert>
#include <stdexcept>
#include <iostream>


namespace robotoc {

LineSearchFilter::LineSearchFilter(const double beta) 
  : filter_(),
    beta_(beta) {
  if (beta_ <= 0) {
    throw std::out_of_range("[LineSearchFilter] invalid argument: beta must be positive!");
  }
}


bool LineSearchFilter::isAccepted(const double cost, 
                                  const double constraint_violation) const {
  assert(constraint_violation >= 0);
  if (filter_.empty()) {
    return true;
  }
  for (const auto& e : filter_) {
    if ((cost < e.first - beta_ * e.second) 
        || (constraint_violation < e.second - beta_ * e.second)) {
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