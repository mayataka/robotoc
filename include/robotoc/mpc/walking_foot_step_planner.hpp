#ifndef ROBOTOC_WALKING_FOOT_STEP_PLANNER_HPP_
#define ROBOTOC_WALKING_FOOT_STEP_PLANNER_HPP_

#include <vector>
#include <iostream>

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/se3.hpp"
#include "robotoc/utils/aligned_vector.hpp"


namespace robotoc {

///
/// @class WalkingFootStepPlanner
/// @brief Foot step planner for the walking gait of biped robot. 
///
class WalkingFootStepPlanner {
public:
  ///
  /// @brief Constructs the planner.
  /// @param[in] biped_robot Biped robot model. 
  ///
  WalkingFootStepPlanner(const Robot& biped_robot);

  ///
  /// @brief Default constructor. 
  ///
  WalkingFootStepPlanner();

  ///
  /// @brief Destructor. 
  ///
  ~WalkingFootStepPlanner();

  ///
  /// @brief Default copy constructor. 
  ///
  WalkingFootStepPlanner(const WalkingFootStepPlanner&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  WalkingFootStepPlanner& operator=(const WalkingFootStepPlanner&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  WalkingFootStepPlanner(WalkingFootStepPlanner&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  WalkingFootStepPlanner& operator=(WalkingFootStepPlanner&&) noexcept = default;

  ///
  /// @brief Sets the gait pattern. 
  /// @param[in] step_length Step length of the gait. 
  /// @param[in] yaw_rate Yaw-rate of the gait. 
  /// @param[in] enable_double_support_phase Falgs to enable the double support
  /// phase or not. 
  ///
  void setGaitPattern(const Eigen::Vector3d& step_length, const double yaw_rate,
                      const bool enable_double_support_phase);

  ///
  /// @brief Initializes the planner. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  ///
  void init(const Eigen::VectorXd& q);

  ///
  /// @brief Plans the foot steps. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  /// @param[in] contact_status Initial contact status.
  /// @param[in] planning_steps Number of planning steps. Must be non-negative.
  /// @return True if the planning is succeeded. False if not.
  ///
  bool plan(const Eigen::VectorXd& q, const ContactStatus& contact_status,
            const int planning_steps);

  ///
  /// @brief Gets the contact placements of a specified step. 
  /// @param[in] step Step of interest.
  /// @return Contact placements of a specified step. 
  ///
  const aligned_vector<SE3>& contactPlacement(const int step) const;

  ///
  /// @brief Gets the CoM position of a specified step. 
  /// @param[in] step Step of interest.
  /// @return CoM position of a specified step. 
  ///
  const Eigen::Vector3d& com(const int step) const;

  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const WalkingFootStepPlanner& planner);

  friend std::ostream& operator<<(std::ostream& os, 
                                  const std::shared_ptr<WalkingFootStepPlanner>& planner);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Robot robot_;
  int L_foot_id_, R_foot_id_, previous_initial_step_;
  double left_to_right_leg_distance_;
  aligned_vector<aligned_vector<SE3>> contact_placement_ref_;
  std::vector<Eigen::Vector3d> com_ref_, com_to_contact_position_local_;
  Eigen::Vector3d step_length_, com_;
  Eigen::Matrix3d R_yaw_;
  bool enable_double_support_phase_;

};

} // namespace robotoc 

#endif // ROBOTOC_WALKING_FOOT_STEP_PLANNER_HPP_ 