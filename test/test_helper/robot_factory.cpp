#include "robot_factory.hpp"
#include "urdf_factory.hpp"

#include <string>

#include "robotoc/robot/robot.hpp"


namespace robotoc {
namespace testhelper {

Robot CreateRobotManipulator() {
  const std::string urdf = RobotManipulatorURDF();
  return Robot(RobotModelInfo::Manipulator(urdf));
}


Robot CreateRobotManipulator(const double time_step) {
  assert(time_step >= 0);
  const std::string urdf = RobotManipulatorURDF();
  auto info = RobotModelInfo::Manipulator(urdf);
  info.point_contacts.push_back(ContactModelInfo("iiwa_link_ee_kuka", time_step));
  return Robot(info);
}


Robot CreateQuadrupedalRobot() {
  const std::string urdf = QuadrupedURDF();
  return Robot(RobotModelInfo::Quadruped(urdf, {}));
}


Robot CreateQuadrupedalRobot(const double time_step) {
  assert(time_step >= 0);
  const std::string urdf = QuadrupedURDF();
  const std::vector<std::string> contact_frames = {"LF_FOOT", "LH_FOOT", "RF_FOOT", "RH_FOOT"}; 
  std::vector<ContactModelInfo> point_contacts;
  for (const auto& e : contact_frames) {
    point_contacts.push_back(ContactModelInfo(e, time_step));
  }
  return Robot(RobotModelInfo::Quadruped(urdf, point_contacts));
}


Robot CreateHumanoidRobot() {
  const std::string urdf = HumanoidURDF();
  return Robot(RobotModelInfo::Humanoid(urdf, {}));
}


Robot CreateHumanoidRobot(const double time_step) {
  assert(time_step >= 0);
  const std::string urdf = HumanoidURDF();
  const std::vector<std::string> contact_frames = {"l_sole", "r_sole"}; 
  std::vector<ContactModelInfo> surface_contacts;
  for (const auto& e : contact_frames) {
    surface_contacts.push_back(ContactModelInfo(e, time_step));
  }
  return Robot(RobotModelInfo::Humanoid(urdf, surface_contacts));
}

} // namespace testhelper
} // namespace robotoc