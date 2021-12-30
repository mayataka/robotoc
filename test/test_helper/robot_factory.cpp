#include "robot_factory.hpp"

#include <string>

#include "robotoc/robot/robot.hpp"


namespace robotoc {
namespace testhelper {

Robot CreateRobotManipulator() {
  const std::string urdf = "../urdf/iiwa14/iiwa14.urdf";
  return Robot(urdf);
}


Robot CreateRobotManipulator(const double time_step) {
  assert(time_step >= 0);
  const std::string urdf = "../urdf/iiwa14/iiwa14.urdf";
  const std::vector<int> contact_frames = {18};
  const std::vector<ContactType> contact_types(contact_frames.size(), ContactType::PointContact);
  return Robot(urdf, BaseJointType::FixedBase, contact_frames, contact_types, time_step);
}


Robot CreateQuadrupedalRobot() {
  const std::string urdf = "../urdf/anymal/anymal.urdf";
  return Robot(urdf, BaseJointType::FloatingBase);
}


Robot CreateQuadrupedalRobot(const double time_step) {
  assert(time_step >= 0);
  const std::string urdf = "../urdf/anymal/anymal.urdf";
  const std::vector<int> contact_frames = {12, 22, 32, 42};
  const std::vector<ContactType> contact_types(contact_frames.size(), ContactType::PointContact);
  return Robot(urdf, BaseJointType::FloatingBase, contact_frames, contact_types, time_step);
}


Robot CreateHumanoidRobot() {
  const std::string urdf = "../urdf/icub/icub.urdf";
  return Robot(urdf, BaseJointType::FloatingBase);
}


Robot CreateHumanoidRobot(const double time_step) {
  assert(time_step >= 0);
  const std::string urdf = "../urdf/icub/icub.urdf";
  const std::vector<int> contact_frames = {26, 46};
  const std::vector<ContactType> contact_types(contact_frames.size(), ContactType::SurfaceContact);
  return Robot(urdf, BaseJointType::FloatingBase, contact_frames, contact_types, time_step);
}

} // namespace testhelper
} // namespace robotoc