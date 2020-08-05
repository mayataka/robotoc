#ifndef IDOCP_COST_FUNCTION_FACTORY_HPP_
#define IDOCP_COST_FUNCTION_FACTORY_HPP_

#include <memory>

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_interface.hpp"
#include "idocp/cost/cost_function_factory_interface.hpp"


namespace idocp {

template <class CostFunctionType>
class CostFunctionFactory final : public CostFunctionFactoryInterface {
public:
  CostFunctionFactory();

  ~CostFunctionFactory();

  // Use default copy constructor.
  CostFunctionFactory(const CostFunctionFactory&) = default;

  // Use default copy coperator.
  CostFunctionFactory& operator=(const CostFunctionFactory&) = default;

  // Use default move constructor.
  CostFunctionFactory(CostFunctionFactory&&) noexcept = default;

  // Use default move assign coperator.
  CostFunctionFactory& operator=(CostFunctionFactory&&) noexcept = default;

  std::unique_ptr<CostFunctionInterface> create(const Robot& robot) override {
    return std::make_unique<CostFunctionType>(robot);
  }
};

} // namespace idocp


#endif // IDOCP_COST_FUNCTION_FACTORY_HPP_