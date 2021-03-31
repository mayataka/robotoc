#ifndef IDOCP_TEST_HELPER_COST_FACTORY_HPP_
#define IDOCP_TEST_HELPER_COST_FACTORY_HPP_

#include <memory>

#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/time_varying_configuration_space_cost.hpp"


namespace idocp {
namespace testhelper {
class TimeVaryingConfigurationRef : public TimeVaryingConfigurationRefBase {
public:
  TimeVaryingConfigurationRef(const Eigen::VectorXd& q0_ref, 
                              const Eigen::VectorXd& v_ref);

  TimeVaryingConfigurationRef() {}

  ~TimeVaryingConfigurationRef() {}

  TimeVaryingConfigurationRef(const TimeVaryingConfigurationRef&) = default;

  TimeVaryingConfigurationRef& operator=( 
      const TimeVaryingConfigurationRef&) = default;

  TimeVaryingConfigurationRef(
      TimeVaryingConfigurationRef&&) noexcept = default;

  TimeVaryingConfigurationRef& operator=(
      TimeVaryingConfigurationRef&&) noexcept = default;

  void update_q_ref(const Robot& robot, const double t, 
                    Eigen::VectorXd& q_ref) const override;

  bool isActive(const double t) const override;

private:
  Eigen::VectorXd q0_ref_, v_ref_;
};


std::shared_ptr<CostFunction> CreateCost(const Robot& robot);

} // namespace testhelper
} // namespace idocp

#endif // IDOCP_TEST_HELPER_COST_FACTORY_HPP_