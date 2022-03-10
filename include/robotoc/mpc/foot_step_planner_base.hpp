#ifndef ROBOTOC_FOOT_STEP_PLANNER_BASE_HPP_
#define ROBOTOC_FOOT_STEP_PLANNER_BASE_HPP_

#include <vector>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/se3.hpp"
#include "robotoc/utils/aligned_vector.hpp"


namespace robotoc {

///
/// @class FootStepPlannerBase
/// @brief Base interface of foot step planners.
///
class FootStepPlannerBase {
public:
  ///
  /// @brief Default constructor. 
  ///
  FootStepPlannerBase() {}

  ///
  /// @brief Destructor. 
  ///
  virtual ~FootStepPlannerBase() {}

  ///
  /// @brief Default copy constructor. 
  ///
  FootStepPlannerBase(const FootStepPlannerBase&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  FootStepPlannerBase& operator=(const FootStepPlannerBase&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  FootStepPlannerBase(FootStepPlannerBase&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  FootStepPlannerBase& operator=(FootStepPlannerBase&&) noexcept = default;

  ///
  /// @brief Initializes the planner. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  ///
  virtual void init(const Eigen::VectorXd& q) = 0;

  ///
  /// @brief Plans the foot steps. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  /// @param[in] contact_status Initial contact status.
  /// @param[in] planning_steps Number of planning steps. Must be non-negative.
  /// @return True if the planning is succeeded. False if not.
  ///
  virtual bool plan(const Eigen::VectorXd& q, 
                    const ContactStatus& contact_status,
                    const int planning_steps) = 0;

  ///
  /// @brief Gets the contact placements of a specified step. 
  /// @param[in] step Step of interest.
  /// @return Contact placements of a specified step. 
  ///
  virtual const aligned_vector<SE3>& contactPlacement(const int step) const = 0;

  ///
  /// @brief Gets the contact positions of a specified step. 
  /// @param[in] step Step of interest.
  /// @return Contact positions of a specified step. 
  ///
  virtual const std::vector<Eigen::Vector3d>& contactPosition(const int step) const = 0;

  ///
  /// @brief Gets the CoM position of a specified step. 
  /// @param[in] step Step of interest.
  /// @return CoM position of a specified step. 
  ///
  virtual const Eigen::Vector3d& com(const int step) const = 0;

};

} // namespace robotoc 

#endif // ROBOTOC_FOOT_STEP_PLANNER_BASE_HPP_ 