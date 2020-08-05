#ifndef IDOCP_MANIPULATOR_CONSTRAINTS_FACTORY_HPP_
#define IDOCP_MANIPULATOR_CONSTRAINTS_FACTORY_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/constraints/constraints_factory_interface.hpp"


namespace idocp {
namespace manipulator {

class ConstraintsFactory final : public ConstraintsFactoryInterface {
public:
  ConstraintsFactory();

  ~ConstraintsFactory();

  ConstraintsFactory(const ConstraintsFactory&) = default;

  ConstraintsFactory& operator=(const ConstraintsFactory&) = default;

  ConstraintsFactory(ConstraintsFactory&&) noexcept = default;

  ConstraintsFactory& operator=(ConstraintsFactory&&) noexcept = default;

  std::unique_ptr<ConstraintsInterface> create(const Robot& robot) override;
};

} // namespace manipulator
} // namespace idocp


#endif // IDOCP_MANIPULATOR_CONSTRAINTS_FACTORY_HPP_ 