#include "robotoc/robot/robot_model_info.hpp"

namespace robotoc {

RobotModelInfo::RobotModelInfo(const std::string& _urdf_path, 
                               const BaseJointType _base_joint_type,
                               const std::vector<ContactModelInfo>& _point_contacts,
                               const std::vector<ContactModelInfo>& _surface_contacts,
                               const double _contact_inv_damping)
  : urdf_path(_urdf_path),
    base_joint_type(_base_joint_type),
    point_contacts(_point_contacts),
    surface_contacts(_surface_contacts),
    contact_inv_damping(_contact_inv_damping) {
}


RobotModelInfo RobotModelInfo::Manipulator(const std::string& urdf_path) {
  return RobotModelInfo(urdf_path, BaseJointType::FixedBase, {}, {});
}


RobotModelInfo RobotModelInfo::Quadruped(
    const std::string& urdf_path, 
    const std::vector<ContactModelInfo>& point_contacts) {
  return RobotModelInfo(urdf_path, BaseJointType::FloatingBase, point_contacts, {});
}


RobotModelInfo RobotModelInfo::Humanoid(
    const std::string& urdf_path, 
    const std::vector<ContactModelInfo>& surface_contacts) {
  return RobotModelInfo(urdf_path, BaseJointType::FloatingBase, {}, surface_contacts);
}

} // namespace robotoc