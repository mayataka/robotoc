#ifndef ROBOTOC_BIPED_WALK_FOOT_STEP_PLANNER_HPP_
#define ROBOTOC_BIPED_WALK_FOOT_STEP_PLANNER_HPP_

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
/// @class BipedWalkFootStepPlanner
/// @brief Foot step planner for the biped robot walk. 
///
class BipedWalkFootStepPlanner final : public ContactPlannerBase {
public:
  ///
  /// @brief Constructs the planner.
  /// @param[in] biped_robot Biped robot model. 
  ///
  BipedWalkFootStepPlanner(const Robot& biped_robot);

  ///
  /// @brief Default constructor. 
  ///
  BipedWalkFootStepPlanner();

  ///
  /// @brief Destructor. 
  ///
  ~BipedWalkFootStepPlanner();

  ///
  /// @brief Default copy constructor. 
  ///
  BipedWalkFootStepPlanner(const BipedWalkFootStepPlanner&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  BipedWalkFootStepPlanner& operator=(const BipedWalkFootStepPlanner&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  BipedWalkFootStepPlanner(BipedWalkFootStepPlanner&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  BipedWalkFootStepPlanner& operator=(BipedWalkFootStepPlanner&&) noexcept = default;

  DEFINE_DEFAULT_CLONE_CONTACT_PLANNER(BipedWalkFootStepPlanner);

  ///
  /// @brief Sets the gait pattern by step length and yaw step command. 
  /// @param[in] step_length Step length of the gait. 
  /// @param[in] step_yaw Yaw command at each step of the gait. 
  /// @param[in] enable_double_support_phase Enables the double-support 
  /// phase or not. 
  ///
  void setGaitPattern(const Eigen::Vector3d& step_length, const double step_yaw,
                      const bool enable_double_support_phase);

  ///
  /// @brief Sets the gait pattern by Raibert heuristic. 
  /// @param[in] vcom_cmd Command of the COM velocity. 
  /// @param[in] yaw_rate_cmd Command of the yaw-rate of the body. 
  /// @param[in] swing_time Swing time of the gait. 
  /// @param[in] double_support_time Double support time of the gait. 
  /// @param[in] gain The feedback gain of the vcom_cmd. 
  ///
  void setRaibertGaitPattern(const Eigen::Vector3d& vcom_cmd, 
                             const double yaw_rate_cmd, const double swing_time, 
                             const double double_support_time, const double gain);

  void init(const Eigen::VectorXd& q) override;

  bool plan(const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
            const ContactStatus& contact_status, 
            const int planning_steps) override;

  const aligned_vector<SE3>& contactPlacements(const int step) const override;

  ///
  /// @brief This is invalid in BipedWalkFootStepPlanner. 
  ///
  const std::vector<Eigen::Vector3d>& contactPositions(const int step) const override;

  ///
  /// @brief This is invalid in BipedWalkFootStepPlanner. 
  ///
  const std::vector<Eigen::Matrix3d>& contactSurfaces(const int step) const override;

  const Eigen::Vector3d& CoM(const int step) const override;

  const Eigen::Matrix3d& R(const int step) const override;

  int size() const override { return planning_size_; }

  friend std::ostream& operator<<(std::ostream& os, 
                                  const BipedWalkFootStepPlanner& planner);

  friend std::ostream& operator<<(std::ostream& os, 
                                  const std::shared_ptr<BipedWalkFootStepPlanner>& planner);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Robot robot_;
  RaibertHeuristic raibert_heuristic_;
  MovingWindowFilter<2> vcom_moving_window_filter_;
  bool enable_raibert_heuristic_;
  int L_foot_id_, R_foot_id_, current_step_, planning_size_;
  double left_to_right_leg_distance_, foot_height_to_com_height_;
  aligned_vector<aligned_vector<SE3>> contact_placement_ref_;
  std::vector<std::vector<Eigen::Vector3d>> contact_position_ref_;
  std::vector<std::vector<Eigen::Matrix3d>> contact_surface_ref_;
  std::vector<Eigen::Vector3d> com_ref_;
  std::vector<Eigen::Matrix3d> R_;
  Eigen::Vector3d vcom_, vcom_cmd_, step_length_;
  Eigen::Matrix3d R_yaw_, R_current_;
  double yaw_rate_cmd_;
  bool enable_double_support_phase_;

};

} // namespace robotoc 

#endif // ROBOTOC_BIPED_WALK_FOOT_STEP_PLANNER_HPP_