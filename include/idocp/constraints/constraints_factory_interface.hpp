#ifndef IDOCP_CONSTRAINTS_FACTORY_INTERFACE_HPP_
#define IDOCP_CONSTRAINTS_FACTORY_INTERFACE_HPP_

#include <memory>

#include "idocp/robot/robot.hpp"
#include "idocp/constraints/constraints_interface.hpp"


namespace idocp {

class ConstraintsFactoryInterface {
public:
  ConstraintsFactoryInterface() {}

  virtual ~ConstraintsFactoryInterface() {}

  // Use default copy constructor.
  ConstraintsFactoryInterface(const ConstraintsFactoryInterface&) = default;

  // Use default copy coperator.
  ConstraintsFactoryInterface& operator=(const ConstraintsFactoryInterface&) = default;

  // Use default move constructor.
  ConstraintsFactoryInterface(ConstraintsFactoryInterface&&) noexcept = default;

  // Use default move assign coperator.
  ConstraintsFactoryInterface& operator=(ConstraintsFactoryInterface&&) noexcept = default;

  virtual std::unique_ptr<ConstraintsInterface> create(const Robot& robot) = 0;
};

} // namespace idocp

#endif // IDOCP_CONSTRAINTS_FACTORY_INTERFACE_HPP_