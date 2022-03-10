#ifndef ROBOTOC_TROTTING_FOOT_STEP_PLANNER_HPP_
#define ROBOTOC_TROTTING_FOOT_STEP_PLANNER_HPP_

#include <vector>
#include <iostream>

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"


namespace robotoc {

///
/// @class TrottingFootStepPlanner
/// @brief Foot step planner for the trotting gait of quadrupeds. 
///
class TrottingFootStepPlanner {
public:
  ///
  /// @brief Constructs the planner.
  /// @param[in] quadruped_robot Quadruped robot model. 
  ///
  TrottingFootStepPlanner(const Robot& quadruped_robot);

  ///
  /// @brief Default constructor. 
  ///
  TrottingFootStepPlanner();

  ///
  /// @brief Destructor. 
  ///
  ~TrottingFootStepPlanner();

  ///
  /// @brief Default copy constructor. 
  ///
  TrottingFootStepPlanner(const TrottingFootStepPlanner&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  TrottingFootStepPlanner& operator=(const TrottingFootStepPlanner&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  TrottingFootStepPlanner(TrottingFootStepPlanner&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  TrottingFootStepPlanner& operator=(TrottingFootStepPlanner&&) noexcept = default;

  ///
  /// @brief Sets the gait pattern. 
  /// @param[in] step_length Step length of the gait. 
  /// @param[in] yaw_rate Yaw-rate of the gait. 
  ///
  void setGaitPattern(const Eigen::Vector3d& step_length, const double yaw_rate);

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
  /// @brief Gets the contact positions of a specified step. 
  /// @param[in] step Step of interest.
  /// @return Contact positions of a specified step. 
  ///
  const std::vector<Eigen::Vector3d>& contactPosition(const int step) const;

  ///
  /// @brief Gets the CoM position of a specified step. 
  /// @param[in] step Step of interest.
  /// @return CoM position of a specified step. 
  ///
  const Eigen::Vector3d& com(const int step) const;

  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const TrottingFootStepPlanner& planner);

  friend std::ostream& operator<<(std::ostream& os, 
                                  const std::shared_ptr<TrottingFootStepPlanner>& planner);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Robot robot_;
  int LF_foot_id_, LH_foot_id_, RF_foot_id_, RH_foot_id_;
  std::vector<std::vector<Eigen::Vector3d>> contact_position_ref_;
  std::vector<Eigen::Vector3d> com_ref_, com_to_contact_position_local_;
  Eigen::Vector3d step_length_, com_;
  Eigen::Matrix3d R_yaw_;

};

} // namespace robotoc 

#endif // ROBOTOC_TROTTING_FOOT_STEP_PLANNER_HPP_ 