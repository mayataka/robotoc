#ifndef ROBOTOC_TEST_HELPER_ROBOT_FACTORY_HPP_
#define ROBOTOC_TEST_HELPER_ROBOT_FACTORY_HPP_

#include "robotoc/robot/robot.hpp"


namespace robotoc {
namespace testhelper {

Robot CreateFixedBaseRobot();

Robot CreateFixedBaseRobot(const double time_step);

Robot CreateFloatingBaseRobot();

Robot CreateFloatingBaseRobot(const double time_step);

} // namespace testhelper
} // namespace robotoc

#endif // ROBOTOC_TEST_HELPER_ROBOT_FACTORY_HPP_ 