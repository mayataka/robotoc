#ifndef IDOCP_PERIODIC_SWITCHING_TIME_COST_HPP_
#define IDOCP_PERIODIC_SWITCHING_TIME_COST_HPP_

#include "idocp/cost/switching_time_cost_function_component_base.hpp"

#include "Eigen/Core"


namespace idocp {

///
/// @class PeriodicSwitchingTimeCost
/// @brief Cost on the deviation from the periodic reference switching time.
///
class PeriodicSwitchingTimeCost : public SwitchingTimeCostFunctionComponentBase {
public:
  PeriodicSwitchingTimeCost(const double period, const double t_start);

  PeriodicSwitchingTimeCost();

  ~PeriodicSwitchingTimeCost();

  PeriodicSwitchingTimeCost(const PeriodicSwitchingTimeCost&) = default;

  PeriodicSwitchingTimeCost& operator=(
      const PeriodicSwitchingTimeCost&) = default;

  PeriodicSwitchingTimeCost(PeriodicSwitchingTimeCost&&) noexcept = default;

  PeriodicSwitchingTimeCost& operator=(
      PeriodicSwitchingTimeCost&&) noexcept = default;

  void set_period(const double period);

  double computeCost(const double t0, const double tf, 
                     const Eigen::VectorXd& ts) const override;

  void computeCostDerivatives(const double t0, const double tf, 
                              const Eigen::VectorXd& ts, 
                              Eigen::VectorXd& hts) const override;

  void computeCostHessian(const double t0, const double tf, 
                          const Eigen::VectorXd& ts,
                          Eigen::MatrixXd& Qts) const override;

private:
  double period_, t_start_;

};

} // namespace idocp

#endif // IDOCP_PERIODIC_SWITCHING_TIME_COST_HPP_ 