#ifndef IDOCP_MANIPULATOR_COST_FUNCTION_FACTORY_HPP_
#define IDOCP_MANIPULATOR_COST_FUNCTION_FACTORY_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_factory_interface.hpp"


namespace idocp {
namespace manipulator {

class CostFunctionFactory final : public CostFunctionFactoryInterface {
public:
  CostFunctionFactory();

  ~CostFunctionFactory();

  CostFunctionFactory(const CostFunctionFactory&) = default;

  CostFunctionFactory& operator=(const CostFunctionFactory&) = default;

  CostFunctionFactory(CostFunctionFactory&&) noexcept = default;

  CostFunctionFactory& operator=(CostFunctionFactory&&) noexcept = default;

  std::unique_ptr<CostFunctionInterface> create(const Robot& robot) override;
};

} // namespace manipulator
} // namespace idocp


#endif // IDOCP_MANIPULATOR_COST_FUNCTION_FACTORY_HPP_