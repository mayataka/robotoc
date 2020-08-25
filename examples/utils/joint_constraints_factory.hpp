#ifndef IDOCP_UTILS_JOINT_CONSTRAINTS_FACTORY_HPP_
#define IDOCP_UTILS_JOINT_CONSTRAINTS_FACTORY_HPP_

#include <memory>

#include "idocp/robot/robot.hpp"

namespace idocp {

class JointConstraintsFactory {
public:
  JointConstraintsFactory(const Robot& robot);

  std::shared_ptr<idocp::Constraints> create() const;

private:
  Robot robot_;

};

} // namespace idocp

#include "joint_constraints_factory.hxx"

#endif // IDOCP_UTILS_JOINT_CONSTRAINTS_FACTORY_HPP_