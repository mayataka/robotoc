#ifndef IDOCP_UTILS_JOINT_CONSTRAINTS_FACTORY_HPP_
#define IDOCP_UTILS_JOINT_CONSTRAINTS_FACTORY_HPP_

#include <memory>

#include "idocp/robot/robot.hpp"
#include "idocp/constraints/constraints.hpp"

namespace idocp {

class JointConstraintsFactory {
public:
  JointConstraintsFactory(const Robot& robot);

  ~JointConstraintsFactory();

  std::shared_ptr<idocp::Constraints> create() const;

private:
  Robot robot_;

};

} // namespace idocp

#endif // IDOCP_UTILS_JOINT_CONSTRAINTS_FACTORY_HPP_