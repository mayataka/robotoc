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
#include "robotoc/mpc/contact_planner_base.hpp"
#include "robotoc/mpc/raibert_heuristic.hpp"


namespace robotoc {

///
/// @class WalkingFootStepPlanner
/// @brief Foot step planner for the walking gait of biped robot. 
///
class WalkingFootStepPlanner final : public ContactPlannerBase {
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
  /// @param[in] v_com_cmd Command of the COM velocity. 
  /// @param[in] yaw_rate_cmd Command of the yaw-rate of the body. 
  /// @param[in] t_swing Duration of swing. 
  /// @param[in] t_stance Duration of stance. 
  /// @param[in] gain The feedback gain of the v_com_cmd. 
  ///
  void setGaitPattern(const Eigen::Vector3d& v_com_cmd, 
                      const double yaw_rate_cmd, const double t_swing, 
                      const double t_stance, const double gain);

  void init(const Eigen::VectorXd& q) override;

  bool plan(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
            const ContactStatus& contact_status, 
            const int planning_steps) override;

  const aligned_vector<SE3>& contactPlacement(const int step) const override;

  const aligned_vector<aligned_vector<SE3>>& contactPlacement() const override;

  const std::vector<Eigen::Vector3d>& contactPosition(const int step) const override;

  const std::vector<std::vector<Eigen::Vector3d>>& contactPosition() const override;

  const Eigen::Vector3d& com(const int step) const override;

  const std::vector<Eigen::Vector3d>& com() const override;

  const Eigen::Matrix3d& R(const int step) const override;

  const std::vector<Eigen::Matrix3d>& R() const override;

  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const WalkingFootStepPlanner& planner);

  friend std::ostream& operator<<(std::ostream& os, 
                                  const std::shared_ptr<WalkingFootStepPlanner>& planner);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Robot robot_;
  RaibertHeuristic raibert_heuristic_;
  bool enable_raibert_heuristic_;
  int L_foot_id_, R_foot_id_, current_step_;
  double left_to_right_leg_distance_, foot_height_to_com_height_;
  aligned_vector<aligned_vector<SE3>> contact_placement_ref_;
  std::vector<std::vector<Eigen::Vector3d>> contact_position_ref_;
  std::vector<Eigen::Vector3d> com_ref_;
  std::vector<Eigen::Matrix3d> R_;
  Eigen::Vector3d v_com_cmd_, step_length_;
  Eigen::Matrix3d R_yaw_;
  double yaw_rate_cmd_;
  bool enable_double_support_phase_;

};

} // namespace robotoc 

#endif // ROBOTOC_WALKING_FOOT_STEP_PLANNER_HPP_ 