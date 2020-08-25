#ifndef IDOCP_UTILS_JOINT_CONSTRAINTS_FACTORY_HXX_
#define IDOCP_UTILS_JOINT_CONSTRAINTS_FACTORY_HXX_

#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"

namespace idocp {

inline JointConstraintsFactory::JointConstraintsFactory(const Robot& robot) 
  : robot_(robot) {
}


inline std::shared_ptr<idocp::Constraints> JointConstraintsFactory::create() const {
  auto constraints = std::make_shared<idocp::Constraints>();
  auto joint_position_lower = std::make_shared<idocp::JointPositionLowerLimit>(robot_);
  auto joint_position_upper = std::make_shared<idocp::JointPositionUpperLimit>(robot_);
  auto joint_velocity_lower = std::make_shared<idocp::JointVelocityLowerLimit>(robot_);
  auto joint_velocity_upper = std::make_shared<idocp::JointVelocityUpperLimit>(robot_);
  auto joint_torques_lower = std::make_shared<idocp::JointTorquesLowerLimit>(robot_);
  auto joint_torques_upper = std::make_shared<idocp::JointTorquesUpperLimit>(robot_);
  constraints->push_back(joint_position_lower);
  constraints->push_back(joint_position_upper);
  constraints->push_back(joint_velocity_lower);
  constraints->push_back(joint_velocity_upper);
  constraints->push_back(joint_torques_lower);
  constraints->push_back(joint_torques_upper);
  return constraints;
}

} // namespace idocp

#endif // IDOCP_UTILS_JOINT_CONSTRAINTS_FACTORY_HXX_