#include "idocp/cost/periodic_switching_time_cost.hpp"

#include <stdexcept>
#include <cassert>

namespace idocp {

PeriodicSwitchingTimeCost::PeriodicSwitchingTimeCost(const double period, 
                                                     const double t_start)
  : SwitchingTimeCostFunctionComponentBase(),
    period_(0),
    t_start_(0) {
}


PeriodicSwitchingTimeCost::PeriodicSwitchingTimeCost()
  : SwitchingTimeCostFunctionComponentBase(),
    period_(0),
    t_start_(0) {
}


PeriodicSwitchingTimeCost::~PeriodicSwitchingTimeCost() {
}


void PeriodicSwitchingTimeCost::set_period(const double period) {

}


double PeriodicSwitchingTimeCost::computeCost(const double t0, const double tf,
                                              const Eigen::VectorXd& ts) const {
  return 0;
}


void PeriodicSwitchingTimeCost::computeCostDerivatives(
    const double t0, const double tf, const Eigen::VectorXd& ts, 
    Eigen::VectorXd& hts) const {
}


void PeriodicSwitchingTimeCost::computeCostHessian(const double t0, 
                                                   const double tf, 
                                                   const Eigen::VectorXd& ts,
                                                   Eigen::MatrixXd& Qts) const {
}

} // namespace idocp