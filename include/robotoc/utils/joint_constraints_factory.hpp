#ifndef ROBOTOC_UTILS_JOINT_CONSTRAINTS_FACTORY_HPP_
#define ROBOTOC_UTILS_JOINT_CONSTRAINTS_FACTORY_HPP_

#include <memory>

#include "robotoc/robot/robot.hpp"
#include "robotoc/constraints/constraints.hpp"

namespace robotoc {

class JointConstraintsFactory {
public:
  JointConstraintsFactory(const Robot& robot);

  ~JointConstraintsFactory();

  std::shared_ptr<robotoc::Constraints> create() const;

private:
  Robot robot_;

};

} // namespace robotoc

#endif // ROBOTOC_UTILS_JOINT_CONSTRAINTS_FACTORY_HPP_