#ifndef ROBOTOC_TROT_FOOT_STEP_PLANNER_HPP_
#define ROBOTOC_TROT_FOOT_STEP_PLANNER_HPP_

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
/// @class TrotFootStepPlanner
/// @brief Foot step planner for the trot gait of quadrupeds. 
///
class TrotFootStepPlanner final : public ContactPlannerBase {
public:
  ///
  /// @brief Constructs the planner.
  /// @param[in] quadruped_robot Quadruped robot model. 
  ///
  TrotFootStepPlanner(const Robot& quadruped_robot);

  ///
  /// @brief Default constructor. 
  ///
  TrotFootStepPlanner();

  ///
  /// @brief Destructor. 
  ///
  ~TrotFootStepPlanner();

  ///
  /// @brief Default copy constructor. 
  ///
  TrotFootStepPlanner(const TrotFootStepPlanner&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  TrotFootStepPlanner& operator=(const TrotFootStepPlanner&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  TrotFootStepPlanner(TrotFootStepPlanner&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  TrotFootStepPlanner& operator=(TrotFootStepPlanner&&) noexcept = default;

  DEFINE_DEFAULT_CLONE_CONTACT_PLANNER(TrotFootStepPlanner)

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

  ///
  /// @brief Sets the rotation of the contact surfaces. 
  /// @param[in] contact_surfaces Rotation of the contact surfaces. 
  ///
  void setContactSurfaces(const std::vector<Eigen::Matrix3d>& contact_surfaces);

  ///
  /// @brief Sets the rotation of the contact surfaces over the mutiple steps. 
  /// @param[in] contact_surfaces Rotation of the contact surfaces. 
  ///
  void setContactSurfaces(
      const std::vector<std::vector<Eigen::Matrix3d>>& contact_surfaces);

  void init(const Eigen::VectorXd& q) override;

  bool plan(const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
            const ContactStatus& contact_status, 
            const int planning_steps) override;

  ///
  /// @brief This is invalid in TrotFootStepPlanner. 
  ///
  const aligned_vector<SE3>& contactPlacements(const int step) const override;

  const std::vector<Eigen::Vector3d>& contactPositions(const int step) const override;

  const std::vector<Eigen::Matrix3d>& contactSurfaces(const int step) const override;

  const Eigen::Vector3d& CoM(const int step) const override;

  const Eigen::Matrix3d& R(const int step) const override;

  int size() const override { return planning_size_; }

  friend std::ostream& operator<<(std::ostream& os, 
                                  const TrotFootStepPlanner& planner);

  friend std::ostream& operator<<(std::ostream& os, 
                                  const std::shared_ptr<TrotFootStepPlanner>& planner);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Robot robot_;
  RaibertHeuristic raibert_heuristic_;
  MovingWindowFilter<2> vcom_moving_window_filter_;
  bool enable_raibert_heuristic_;
  int LF_foot_id_, LH_foot_id_, RF_foot_id_, RH_foot_id_, current_step_, planning_size_;
  aligned_vector<aligned_vector<SE3>> contact_placement_ref_;
  std::vector<std::vector<Eigen::Vector3d>> contact_position_ref_;
  std::vector<std::vector<Eigen::Matrix3d>> contact_surface_ref_;
  std::vector<Eigen::Vector3d> com_ref_, com_to_contact_position_local_;
  std::vector<Eigen::Matrix3d> R_;
  Eigen::Vector3d vcom_, vcom_cmd_, step_length_;
  Eigen::Matrix3d R_yaw_, R_current_;
  double yaw_rate_cmd_;
  bool enable_stance_phase_;

};

} // namespace robotoc 

#endif // ROBOTOC_TROT_FOOT_STEP_PLANNER_HPP_ 