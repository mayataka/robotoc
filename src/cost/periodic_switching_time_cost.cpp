#include "idocp/cost/periodic_switching_time_cost.hpp"

#include <stdexcept>
#include <cassert>
#include <iostream>

namespace idocp {

PeriodicSwitchingTimeCost::PeriodicSwitchingTimeCost(const double period, 
                                                     const double t_start)
  : SwitchingTimeCostFunctionComponentBase(),
    period_(0),
    t_start_(0),
    weight_(0) {
  try {
    if (period <= 0) {
      throw std::out_of_range(
          "Invalid argment: period must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


PeriodicSwitchingTimeCost::PeriodicSwitchingTimeCost()
  : SwitchingTimeCostFunctionComponentBase(),
    period_(0),
    t_start_(0) {
}


PeriodicSwitchingTimeCost::~PeriodicSwitchingTimeCost() {
}


void PeriodicSwitchingTimeCost::set_period(const double period,
                                           const double t_start) {
  try {
    if (period <= 0) {
      throw std::out_of_range(
          "Invalid argment: period must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  period_ = period;
  t_start_ = t_start;
}


void PeriodicSwitchingTimeCost::set_weight(const double weight) {
  try {
    if (weight < 0) {
      throw std::out_of_range(
          "Invalid argment: weight must be non-negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  weight_ = weight;
}


double PeriodicSwitchingTimeCost::computeCost(const double t0, const double tf,
                                              const Eigen::VectorXd& ts) const {
  const int num_switches = ts.size();
  double cost = 0;
  for (int i=0; num_switches; ++i) {
    const double ts_ref = t_start_ + (i+1)*period_;
    cost += 0.5 * weight_ * (ts.coeff(i)-ts_ref) * (ts.coeff(i)-ts_ref);
  }
  return cost;
}


void PeriodicSwitchingTimeCost::computeCostDerivatives(
    const double t0, const double tf, const Eigen::VectorXd& ts, 
    Eigen::VectorXd& hts) const {
  const int num_switches = ts.size();
  for (int i=0; num_switches; ++i) {
    const double ts_ref = t_start_ + (i+1)*period_;
    hts.coeffRef(i) += weight_ * (ts.coeff(i)-ts_ref);
  }
}


void PeriodicSwitchingTimeCost::computeCostHessian(const double t0, 
                                                   const double tf, 
                                                   const Eigen::VectorXd& ts,
                                                   Eigen::MatrixXd& Qts) const {
  const int num_switches = ts.size();
  for (int i=0; num_switches; ++i) {
    const double ts_ref = t_start_ + (i+1)*period_;
    Qts.coeffRef(i, i) += weight_;
  }
}

} // namespace idocp