#ifndef ROBOTOC_ROBOT_PROPERTIES_HPP_
#define ROBOTOC_ROBOT_PROPERTIES_HPP_

#include <string>
#include <vector>
#include <utility>
#include <iostream>

#include "Eigen/Core"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/container/aligned-vector.hpp"
#include "pinocchio/spatial/force.hpp"

#include "robotoc/robot/se3.hpp"
#include "robotoc/robot/point_contact.hpp"
#include "robotoc/robot/surface_contact.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/utils/aligned_vector.hpp"


namespace robotoc {

///
/// @class RobotProperties
/// @brief Collection of the robot properties, which can change after 
/// constructing robot models.
///
struct RobotProperties {
  Eigen::VectorXd generalized_momentum_bias; 
};

} // namespace robotoc 

#endif // ROBOTOC_ROBOT_PROPERTIES_HPP_