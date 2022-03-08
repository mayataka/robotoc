#ifndef ROBOTOC_CRAWLING_FOOT_STEP_PLANNER_HPP_
#define ROBOTOC_CRAWLING_FOOT_STEP_PLANNER_HPP_

#include <vector>
#include <iostream>

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"


namespace robotoc {

///
/// @class CrawlingFootStepPlanner
/// @brief MPC solver for the crawling gait of quadrupeds. 
///
class CrawlingFootStepPlanner {
public:
  ///
  /// @brief Constructs the planner.
  /// @param[in] robot robot model. 
  ///
  CrawlingFootStepPlanner(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  CrawlingFootStepPlanner();

  ///
  /// @brief Destructor. 
  ///
  ~CrawlingFootStepPlanner();

  ///
  /// @brief Default copy constructor. 
  ///
  CrawlingFootStepPlanner(const CrawlingFootStepPlanner&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  CrawlingFootStepPlanner& operator=(const CrawlingFootStepPlanner&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  CrawlingFootStepPlanner(CrawlingFootStepPlanner&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  CrawlingFootStepPlanner& operator=(CrawlingFootStepPlanner&&) noexcept = default;

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
  /// @brief Gets the contact positions of a specified phase. 
  /// @param[in] phase Phase of interest.
  /// @return Contact positions of a specified phase. 
  ///
  const std::vector<Eigen::Vector3d>& contactPosition(const int phase) const;

  ///
  /// @brief Gets the CoM position of a specified phase. 
  /// @param[in] phase Phase of interest.
  /// @return CoM position of a specified phase. 
  ///
  const Eigen::Vector3d& com(const int phase) const;

  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const CrawlingFootStepPlanner& planner);

  friend std::ostream& operator<<(std::ostream& os, 
                                  const std::shared_ptr<CrawlingFootStepPlanner>& planner);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Robot robot_;
  int LF_foot_id_, LH_foot_id_, RF_foot_id_, RH_foot_id_;
  std::vector<std::vector<Eigen::Vector3d>> contact_position_ref_;
  std::vector<Eigen::Vector3d> com_ref_, com_to_contact_position_local_;
  Eigen::Vector3d step_length_, com_;
  Eigen::Matrix3d R_yaw_;
  bool first_step_;

};

} // namespace robotoc 

#endif // ROBOTOC_CRAWLING_FOOT_STEP_PLANNER_HPP_ 