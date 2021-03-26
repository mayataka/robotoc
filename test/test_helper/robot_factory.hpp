#ifndef IDOCP_TEST_HELPER_ROBOT_FACTORY_HPP_
#define IDOCP_TEST_HELPER_ROBOT_FACTORY_HPP_

#include "idocp/robot/robot.hpp"


namespace idocp {
namespace testhelper {

Robot CreateFixedBaseRobot();

Robot CreateFixedBaseRobot(const double time_step);

Robot CreateFloatingBaseRobot();

Robot CreateFloatingBaseRobot(const double time_step);

} // namespace testhelper
} // namespace idocp

#endif // IDOCP_TEST_HELPER_ROBOT_FACTORY_HPP_ 