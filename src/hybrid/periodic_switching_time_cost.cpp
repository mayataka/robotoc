#include "robotoc/hybrid/periodic_switching_time_cost.hpp"

#include <stdexcept>
#include <cassert>
#include <iostream>

namespace robotoc {

PeriodicSwitchingTimeCost::PeriodicSwitchingTimeCost(const double period, 
                                                     const double t_start)
  : STOCostFunctionComponentBase(),
    period_(0),
    t_start_(0),
    weight_(0) {
  if (period <= 0) {
    throw std::out_of_range(
        "Invalid argment: period must be positive!");
  }
}


PeriodicSwitchingTimeCost::PeriodicSwitchingTimeCost()
  : STOCostFunctionComponentBase(),
    period_(0),
    t_start_(0) {
}


PeriodicSwitchingTimeCost::~PeriodicSwitchingTimeCost() {
}


void PeriodicSwitchingTimeCost::set_period(const double period,
                                           const double t_start) {
  if (period <= 0) {
    throw std::out_of_range(
        "Invalid argment: period must be positive!");
  }
  period_ = period;
  t_start_ = t_start;
}


void PeriodicSwitchingTimeCost::set_weight(const double weight) {
  if (weight < 0) {
    throw std::out_of_range(
        "Invalid argment: weight must be non-negative!");
  }
  weight_ = weight;
}


double PeriodicSwitchingTimeCost::evalCost(
    const TimeDiscretization& discretization) const {
  double cost = 0;
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  if (num_events > 0) {
    int impulse_index = 0;
    int lift_index = 0;
    if (discretization.eventType(0) == DiscreteEventType::Impulse) {
      const double ts_diff = discretization.impulseTime(0) - t_start_;
      cost += 0.5 * weight_ * ts_diff * ts_diff;
      ++impulse_index;
    }
    else {
      const double ts_diff = discretization.liftTime(0) - t_start_;
      cost += 0.5 * weight_ * ts_diff * ts_diff;
      ++lift_index;
    }
    for (int event_index=1; event_index<num_events; ++event_index) {
      assert(event_index == impulse_index+lift_index);
      if (discretization.eventType(event_index) == DiscreteEventType::Impulse) {
        const double ts_diff = discretization.impulseTime(impulse_index) - t_start_;
        cost += 0.5 * weight_ * ts_diff * ts_diff;
        ++impulse_index;
      }
      else {
        const double ts_diff = discretization.liftTime(lift_index) - t_start_;
        cost += 0.5 * weight_ * ts_diff * ts_diff;
        ++lift_index;
      }
    }
  }
  return cost;
}


void PeriodicSwitchingTimeCost::evalCostDerivatives(
    const TimeDiscretization& discretization, Eigen::VectorXd& lts) const {
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  if (num_events > 0) {
    int impulse_index = 0;
    int lift_index = 0;
    if (discretization.eventType(0) == DiscreteEventType::Impulse) {
      const double ts_diff = discretization.impulseTime(0) - t_start_;
      lts.coeffRef(0) += weight_ * ts_diff;
      ++impulse_index;
    }
    else {
      const double ts_diff = discretization.liftTime(0) - t_start_;
      lts.coeffRef(0) += weight_ * ts_diff;
      ++lift_index;
    }
    for (int event_index=1; event_index<num_events; ++event_index) {
      assert(event_index == impulse_index+lift_index);
      if (discretization.eventType(event_index) == DiscreteEventType::Impulse) {
        const double ts_diff = discretization.impulseTime(impulse_index) - t_start_;
        lts.coeffRef(event_index) += weight_ * ts_diff;
        ++impulse_index;
      }
      else {
        const double ts_diff = discretization.liftTime(lift_index) - t_start_;
        lts.coeffRef(event_index) += weight_ * ts_diff;
        ++lift_index;
      }
    }
  }
}


void PeriodicSwitchingTimeCost::evalCostHessian(
    const TimeDiscretization& discretization, Eigen::MatrixXd& Qts) const {
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  for (int i=0; num_events; ++i) {
    Qts.coeffRef(i, i) += weight_;
  }
}

} // namespace robotoc