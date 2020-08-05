#ifndef IDOCP_COST_FUNCTION_FACTORY_INTERFACE_HPP_
#define IDOCP_COST_FUNCTION_FACTORY_INTERFACE_HPP_

#include <memory>

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_interface.hpp"


namespace idocp {

class CostFunctionFactoryInterface {
public:
  CostFunctionFactoryInterface() {}

  virtual ~CostFunctionFactoryInterface() {}

  // Use default copy constructor.
  CostFunctionFactoryInterface(const CostFunctionFactoryInterface&) = default;

  // Use default copy coperator.
  CostFunctionFactoryInterface& operator=(const CostFunctionFactoryInterface&) = default;

  // Use default move constructor.
  CostFunctionFactoryInterface(CostFunctionFactoryInterface&&) noexcept = default;

  // Use default move assign coperator.
  CostFunctionFactoryInterface& operator=(CostFunctionFactoryInterface&&) noexcept = default;

  virtual std::unique_ptr<CostFunctionInterface> create(const Robot& robot) = 0;
};

} // namespace idocp

#endif // IDOCP_COST_FUNCTION_FACTORY_INTERFACE_HPP_