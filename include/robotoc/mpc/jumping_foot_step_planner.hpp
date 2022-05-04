#ifndef ROBOTOC_JUMPING_FOOT_STEP_PLANNER_HPP_
#define ROBOTOC_JUMPING_FOOT_STEP_PLANNER_HPP_

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
/// @class JumpingFootStepPlanner
/// @brief Foot step planner for the jumping motion. 
///
class JumpingFootStepPlanner final : public ContactPlannerBase {
public:
  ///
  /// @brief Constructs the planner.
  /// @param[in] robot Robot model. 
  ///
  JumpingFootStepPlanner(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  JumpingFootStepPlanner();

  ///
  /// @brief Destructor. 
  ///
  ~JumpingFootStepPlanner();

  ///
  /// @brief Default copy constructor. 
  ///
  JumpingFootStepPlanner(const JumpingFootStepPlanner&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  JumpingFootStepPlanner& operator=(const JumpingFootStepPlanner&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  JumpingFootStepPlanner(JumpingFootStepPlanner&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  JumpingFootStepPlanner& operator=(JumpingFootStepPlanner&&) noexcept = default;

  ///
  /// @brief Sets the gait pattern. 
  /// @param[in] jump_length Jump length. 
  /// @param[in] jump_yaw Change in the yaw angle of the jump. 
  ///
  void setJumpPattern(const Eigen::Vector3d& jump_length, const double jump_yaw);

  void init(const Eigen::VectorXd& q) override;

  bool plan(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
            const ContactStatus& contact_status, 
            const int planning_steps) override;

  const aligned_vector<SE3>& contactPlacement(const int step) const override;

  const aligned_vector<aligned_vector<SE3>>& contactPlacement() const override;

  ///
  /// @brief This is invalid in JumpingFootStepPlanner. 
  ///
  const std::vector<Eigen::Vector3d>& contactPosition(const int step) const override;

  ///
  /// @brief This is invalid in JumpingFootStepPlanner. 
  ///
  const std::vector<std::vector<Eigen::Vector3d>>& contactPosition() const override;

  const Eigen::Vector3d& com(const int step) const override;

  const std::vector<Eigen::Vector3d>& com() const override;

  const Eigen::Matrix3d& R(const int step) const override;

  const std::vector<Eigen::Matrix3d>& R() const override;

  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const JumpingFootStepPlanner& planner);

  friend std::ostream& operator<<(std::ostream& os, 
                                  const std::shared_ptr<JumpingFootStepPlanner>& planner);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Robot robot_;
  std::vector<int> contact_frames_;
  int current_step_;
  aligned_vector<aligned_vector<SE3>> contact_placement_ref_;
  std::vector<std::vector<Eigen::Vector3d>> contact_position_ref_;
  std::vector<Eigen::Vector3d> com_ref_, com_to_contact_position_local_;
  std::vector<Eigen::Matrix3d> R_;
  Eigen::Vector3d jump_length_;
  Eigen::Matrix3d R_yaw_;
  bool is_biped_;

};

} // namespace robotoc 

#endif // ROBOTOC_JUMPING_FOOT_STEP_PLANNER_HPP_ 