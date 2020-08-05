#ifndef IDOCP_CONSTRAINTS_FACTORY_HPP_
#define IDOCP_CONSTRAINTS_FACTORY_HPP_

#include <memory>

#include "idocp/robot/robot.hpp"
#include "idocp/constraints/constraints_interface.hpp"
#include "idocp/constraints/constraints_factory_interface.hpp"


namespace idocp {

template <class ConstraintsType>
class ConstraintsFactory final : public ConstraintsFactoryInterface {
public:
  ConstraintsFactory();

  ~ConstraintsFactory();

  // Use default copy constructor.
  ConstraintsFactory(const ConstraintsFactory&) = default;

  // Use default copy coperator.
  ConstraintsFactory& operator=(const ConstraintsFactory&) = default;

  // Use default move constructor.
  ConstraintsFactory(ConstraintsFactory&&) noexcept = default;

  // Use default move assign coperator.
  ConstraintsFactory& operator=(ConstraintsFactory&&) noexcept = default;

  std::unique_ptr<ConstraintsInterface> create(const Robot& robot) override {
    return std::make_unique<ConstraintsType>(robot);
  }
};

} // namespace idocp


#endif // IDOCP_CONSTRAINTS_FACTORY_HPP_