#include "robotoc/utils/joint_constraints_factory.hpp"

#include "robotoc/constraints/joint_position_lower_limit.hpp"
#include "robotoc/constraints/joint_position_upper_limit.hpp"
#include "robotoc/constraints/joint_velocity_lower_limit.hpp"
#include "robotoc/constraints/joint_velocity_upper_limit.hpp"
#include "robotoc/constraints/joint_torques_lower_limit.hpp"
#include "robotoc/constraints/joint_torques_upper_limit.hpp"


namespace robotoc {

JointConstraintsFactory::JointConstraintsFactory(const Robot& robot) 
  : robot_(robot) {
}


JointConstraintsFactory::~JointConstraintsFactory() {
}


std::shared_ptr<robotoc::Constraints> JointConstraintsFactory::create() const {
  auto constraints = std::make_shared<robotoc::Constraints>();
  auto joint_position_lower = std::make_shared<robotoc::JointPositionLowerLimit>(robot_);
  auto joint_position_upper = std::make_shared<robotoc::JointPositionUpperLimit>(robot_);
  auto joint_velocity_lower = std::make_shared<robotoc::JointVelocityLowerLimit>(robot_);
  auto joint_velocity_upper = std::make_shared<robotoc::JointVelocityUpperLimit>(robot_);
  auto joint_torques_lower = std::make_shared<robotoc::JointTorquesLowerLimit>(robot_);
  auto joint_torques_upper = std::make_shared<robotoc::JointTorquesUpperLimit>(robot_);
  constraints->push_back(joint_position_lower);
  constraints->push_back(joint_position_upper);
  constraints->push_back(joint_velocity_lower);
  constraints->push_back(joint_velocity_upper);
  constraints->push_back(joint_torques_lower);
  constraints->push_back(joint_torques_upper);
  return constraints;
}

} // namespace robotoc