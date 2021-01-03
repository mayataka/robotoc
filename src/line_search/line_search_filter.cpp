#include "idocp/line_search/line_search_filter.hpp"

#include <cassert>
#include <stdexcept>
#include <iostream>


namespace idocp {

LineSearchFilter::LineSearchFilter(const double cost_reduction_rate, 
                                   const double constraints_reduction_rate) 
  : filter_(),
    cost_reduction_rate_(cost_reduction_rate),
    constraints_reduction_rate_(constraints_reduction_rate) {
  try {
    if (cost_reduction_rate_ <= 0) {
      throw std::out_of_range("invalid value: cost_reduction_rate must be positive!");
    }
    if (constraints_reduction_rate_ <= 0) {
      throw std::out_of_range("invalid value: constraints_reduction_rate must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


LineSearchFilter::~LineSearchFilter() {
}


bool LineSearchFilter::isAccepted(const double cost, 
                                  const double constraint_violation) {
  assert(constraint_violation >= 0);
  if (!filter_.empty()) {
    for (auto pair : filter_) {
      if (cost >= pair.first && constraint_violation >= pair.second) {
        return false;
      }
    }
  }
  return true;
}


void LineSearchFilter::augment(const double cost, 
                               const double constraint_violation) {
  assert(constraint_violation >= 0);
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


bool LineSearchFilter::isEmpty() const {
  if (filter_.empty()) {
    return true;
  }
  else {
    return false;
  }
}

} // namespace idocp