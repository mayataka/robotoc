#ifndef ROBOTOC_TEST_HELPER_ROBOT_FACTORY_HPP_
#define ROBOTOC_TEST_HELPER_ROBOT_FACTORY_HPP_

#include "robotoc/robot/robot.hpp"


namespace robotoc {
namespace testhelper {

Robot CreateRobotManipulator();

Robot CreateRobotManipulator(const double time_step);

Robot CreateQuadrupedalRobot();

Robot CreateQuadrupedalRobot(const double time_step);

Robot CreateHumanoidRobot();

Robot CreateHumanoidRobot(const double time_step);

} // namespace testhelper
} // namespace robotoc

#endif // ROBOTOC_TEST_HELPER_ROBOT_FACTORY_HPP_ 