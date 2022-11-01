#include "constraints_factory.hpp"

#include "robotoc/constraints/joint_position_lower_limit.hpp"
#include "robotoc/constraints/joint_position_upper_limit.hpp"
#include "robotoc/constraints/joint_velocity_lower_limit.hpp"
#include "robotoc/constraints/joint_velocity_upper_limit.hpp"
#include "robotoc/constraints/joint_torques_lower_limit.hpp"
#include "robotoc/constraints/joint_torques_upper_limit.hpp"
#include "robotoc/constraints/friction_cone.hpp"
#include "robotoc/constraints/impact_friction_cone.hpp"


namespace robotoc {
namespace testhelper {

std::shared_ptr<Constraints> CreateConstraints(const Robot& robot) {
  auto joint_lower_limit = std::make_shared<JointPositionLowerLimit>(robot);
  auto joint_upper_limit = std::make_shared<JointPositionUpperLimit>(robot);
  auto velocity_lower_limit = std::make_shared<JointVelocityLowerLimit>(robot);
  auto velocity_upper_limit = std::make_shared<JointVelocityUpperLimit>(robot);
  auto torques_lower_limit = std::make_shared<JointTorquesLowerLimit>(robot);
  auto torques_upper_limit = std::make_shared<JointTorquesUpperLimit>(robot);
  auto constraints = std::make_shared<Constraints>();
  constraints->add("joint_upper_limit", joint_upper_limit); 
  constraints->add("joint_lower_limit", joint_lower_limit);
  constraints->add("velocity_lower_limit", velocity_lower_limit); 
  constraints->add("velocity_upper_limit", velocity_upper_limit);
  constraints->add("torques_lower_limit", torques_lower_limit);
  constraints->add("torques_upper_limit", torques_upper_limit);
  if (robot.maxNumContacts() > 0) {
    auto friction_cone = std::make_shared<FrictionCone>(robot);
    auto impact_friction_cone = std::make_shared<ImpactFrictionCone>(robot);
    constraints->add("friction_cone", friction_cone);
    constraints->add("impact_friction_cone", impact_friction_cone);
  }
  return constraints;
}

} // namespace testhelper
} // namespace robotoc