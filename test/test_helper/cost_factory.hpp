#ifndef ROBOTOC_TEST_HELPER_COST_FACTORY_HPP_
#define ROBOTOC_TEST_HELPER_COST_FACTORY_HPP_

#include <memory>

#include "robotoc/cost/cost_function.hpp"
#include "robotoc/cost/configuration_space_ref_base.hpp"


namespace robotoc {
namespace testhelper {
class ConfigurationSpaceRef : public ConfigurationSpaceRefBase {
public:
  ConfigurationSpaceRef(const Eigen::VectorXd& q0_ref, 
                        const Eigen::VectorXd& v_ref);

  ConfigurationSpaceRef() {}

  ~ConfigurationSpaceRef() {}

  ConfigurationSpaceRef(const ConfigurationSpaceRef&) = default;

  ConfigurationSpaceRef& operator=(const ConfigurationSpaceRef&) = default;

  ConfigurationSpaceRef(ConfigurationSpaceRef&&) noexcept = default;

  ConfigurationSpaceRef& operator=(ConfigurationSpaceRef&&) noexcept = default;

  DEFINE_DEFAULT_CLONE_CONFIGURATION_SPACE_REF(ConfigurationSpaceRef)

  void updateRef(const Robot& robot, const GridInfo& grid_info, 
                    Eigen::VectorXd& q_ref) const override;

  bool isActive(const GridInfo& grid_info) const override;

private:
  Eigen::VectorXd q0_ref_, v_ref_;
};


std::shared_ptr<CostFunction> CreateCost(const Robot& robot);

} // namespace testhelper
} // namespace robotoc

#endif // ROBOTOC_TEST_HELPER_COST_FACTORY_HPP_