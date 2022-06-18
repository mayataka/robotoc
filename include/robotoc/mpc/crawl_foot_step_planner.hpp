#ifndef ROBOTOC_CRAWL_FOOT_STEP_PLANNER_HPP_
#define ROBOTOC_CRAWL_FOOT_STEP_PLANNER_HPP_

#include <vector>
#include <iostream>

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/se3.hpp"
#include "robotoc/utils/aligned_vector.hpp"
#include "robotoc/mpc/contact_planner_base.hpp"
#include "robotoc/mpc/raibert_heuristic.hpp"
#include "robotoc/mpc/moving_window_filter.hpp"


namespace robotoc {

///
/// @class CrawlFootStepPlanner
/// @brief Foot step planner for the crawl gait of quadrupeds. 
///
class CrawlFootStepPlanner final : public ContactPlannerBase {
public:
  ///
  /// @brief Constructs the planner.
  /// @param[in] quadruped_robot Quadruped robot model. 
  ///
  CrawlFootStepPlanner(const Robot& quadruped_robot);

  ///
  /// @brief Default constructor. 
  ///
  CrawlFootStepPlanner();

  ///
  /// @brief Destructor. 
  ///
  ~CrawlFootStepPlanner();

  ///
  /// @brief Default copy constructor. 
  ///
  CrawlFootStepPlanner(const CrawlFootStepPlanner&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  CrawlFootStepPlanner& operator=(const CrawlFootStepPlanner&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  CrawlFootStepPlanner(CrawlFootStepPlanner&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  CrawlFootStepPlanner& operator=(CrawlFootStepPlanner&&) noexcept = default;

  ///
  /// @brief Sets the gait pattern by step length and yaw step command. 
  /// @param[in] step_length Step length of the gait. 
  /// @param[in] step_yaw Yaw command at each step of the gait. 
  /// @param[in] enable_stance_phase Enables the stance phase or not. 
  ///
  void setGaitPattern(const Eigen::Vector3d& step_length, const double step_yaw,
                      const bool enable_stance_phase);

  ///
  /// @brief Sets the gait pattern by Raibert heuristic. 
  /// @param[in] vcom_cmd Command of the COM velocity. 
  /// @param[in] yaw_rate_cmd Command of the yaw-rate of the body. 
  /// @param[in] swing_time Swing time of the gait. 
  /// @param[in] stance_time Stance time of the gait. 
  /// @param[in] gain The feedback gain of the vcom_cmd. 
  ///
  void setRaibertGaitPattern(const Eigen::Vector3d& vcom_cmd, 
                             const double yaw_rate_cmd, const double swing_time,
                             const double stance_time, const double gain);

  void init(const Eigen::VectorXd& q) override;

  bool plan(const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
            const ContactStatus& contact_status, 
            const int planning_steps) override;

  ///
  /// @brief This is invalid in CrawlFootStepPlanner. 
  ///
  const aligned_vector<SE3>& contactPlacements(const int step) const override;

  ///
  /// @brief This is invalid in CrawlFootStepPlanner. 
  ///
  const aligned_vector<aligned_vector<SE3>>& contactPlacements() const override;

  const std::vector<Eigen::Vector3d>& contactPositions(const int step) const override;

  const std::vector<std::vector<Eigen::Vector3d>>& contactPositions() const override;

  const Eigen::Vector3d& CoM(const int step) const override;

  const std::vector<Eigen::Vector3d>& CoM() const override;

  const Eigen::Matrix3d& R(const int step) const override;

  const std::vector<Eigen::Matrix3d>& R() const override;

  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const CrawlFootStepPlanner& planner);

  friend std::ostream& operator<<(std::ostream& os, 
                                  const std::shared_ptr<CrawlFootStepPlanner>& planner);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Robot robot_;
  RaibertHeuristic raibert_heuristic_;
  MovingWindowFilter<2> vcom_moving_window_filter_;
  bool enable_raibert_heuristic_;
  int LF_foot_id_, LH_foot_id_, RF_foot_id_, RH_foot_id_, current_step_;
  aligned_vector<aligned_vector<SE3>> contact_placement_ref_;
  std::vector<std::vector<Eigen::Vector3d>> contact_position_ref_;
  std::vector<Eigen::Vector3d> com_ref_, com_to_contact_position_local_;
  std::vector<Eigen::Matrix3d> R_;
  Eigen::Vector3d vcom_, vcom_cmd_, step_length_;
  Eigen::Matrix3d R_yaw_, R_current_;
  double yaw_rate_cmd_;
  bool enable_stance_phase_;
};

} // namespace robotoc 

#endif // ROBOTOC_CRAWL_FOOT_STEP_PLANNER_HPP_ 