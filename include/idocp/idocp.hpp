#ifndef IDOCP_IDOCP_HPP_
#define IDOCP_IDOCP_HPP_

#include "idocp/ocp/parnmpc.hpp"
#include "idocp/ocp/ocp.hpp"
#include "idocp/ocp/mpc.hpp"
#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/contact_force_cost.hpp"
#include "idocp/cost/task_space_3d_cost.hpp"
#include "idocp/cost/task_space_6d_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/constraints_data.hpp"
#include "idocp/constraints/joint_acceleration_lower_limit.hpp"
#include "idocp/constraints/joint_acceleration_upper_limit.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"

#endif // IDOCP_IDOCP_HPP_