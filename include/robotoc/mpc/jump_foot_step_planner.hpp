#ifndef ROBOTOC_JUMP_FOOT_STEP_PLANNER_HPP_
#define ROBOTOC_JUMP_FOOT_STEP_PLANNER_HPP_

#include <vector>
#include <iostream>

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/se3.hpp"
#include "robotoc/utils/aligned_vector.hpp"
#include "robotoc/mpc/contact_planner_base.hpp"


namespace robotoc {

///
/// @class JumpFootStepPlanner
/// @brief Foot step planner for the jump motion. 
///
class JumpFootStepPlanner final : public ContactPlannerBase {
public:
  ///
  /// @brief Constructs the planner.
  /// @param[in] robot Robot model. 
  ///
  JumpFootStepPlanner(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  JumpFootStepPlanner();

  ///
  /// @brief Destructor. 
  ///
  ~JumpFootStepPlanner();

  ///
  /// @brief Default copy constructor. 
  ///
  JumpFootStepPlanner(const JumpFootStepPlanner&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  JumpFootStepPlanner& operator=(const JumpFootStepPlanner&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  JumpFootStepPlanner(JumpFootStepPlanner&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  JumpFootStepPlanner& operator=(JumpFootStepPlanner&&) noexcept = default;

  ///
  /// @brief Sets the gait pattern. 
  /// @param[in] jump_length Jump length. 
  /// @param[in] jump_yaw Change in the yaw angle of the jump. 
  ///
  void setJumpPattern(const Eigen::Vector3d& jump_length, const double jump_yaw);

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

  const aligned_vector<SE3>& contactPlacements(const int step) const override;

  const aligned_vector<aligned_vector<SE3>>& contactPlacements() const override;

  ///
  /// @brief This is invalid in JumpFootStepPlanner. 
  ///
  const std::vector<Eigen::Vector3d>& contactPositions(const int step) const override;

  ///
  /// @brief This is invalid in JumpFootStepPlanner. 
  ///
  const std::vector<std::vector<Eigen::Vector3d>>& contactPositions() const override;

  ///
  /// @brief This is invalid in JumpFootStepPlanner. 
  ///
  const std::vector<Eigen::Matrix3d>& contactSurfaces(const int step) const override;

  ///
  /// @brief This is invalid in JumpFootStepPlanner. 
  ///
  const std::vector<std::vector<Eigen::Matrix3d>>& contactSurfaces() const override;

  const Eigen::Vector3d& CoM(const int step) const override;

  const std::vector<Eigen::Vector3d>& CoM() const override;

  const Eigen::Matrix3d& R(const int step) const override;

  const std::vector<Eigen::Matrix3d>& R() const override;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const JumpFootStepPlanner& planner);

  friend std::ostream& operator<<(std::ostream& os, 
                                  const std::shared_ptr<JumpFootStepPlanner>& planner);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Robot robot_;
  std::vector<int> contact_frames_;
  int current_step_;
  aligned_vector<aligned_vector<SE3>> contact_placement_ref_;
  std::vector<std::vector<Eigen::Vector3d>> contact_position_ref_;
  std::vector<std::vector<Eigen::Matrix3d>> contact_surface_ref_;
  std::vector<Eigen::Vector3d> com_ref_, com_to_contact_position_local_;
  std::vector<Eigen::Matrix3d> R_;
  Eigen::Vector3d jump_length_;
  Eigen::Matrix3d R_yaw_;
  bool is_biped_;

};

} // namespace robotoc 

#endif // ROBOTOC_JUMP_FOOT_STEP_PLANNER_HPP_