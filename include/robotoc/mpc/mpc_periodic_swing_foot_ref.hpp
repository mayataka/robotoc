#ifndef ROBOTOC_MPC_PERIODIC_SWING_FOOT_REF_HPP_
#define ROBOTOC_MPC_PERIODIC_SWING_FOOT_REF_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/task_space_3d_cost.hpp"
#include "robotoc/hybrid/contact_sequence.hpp"
#include "robotoc/mpc/contact_planner_base.hpp"


namespace robotoc {

///
/// @class MPCPeriodicSwingFootRef
/// @brief Periodic reference of the foot position. 
///
class MPCPeriodicSwingFootRef final : public TaskSpace3DRefBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] contact_index Contact index of the foot.
  /// @param[in] swing_height The swing height of the gait.
  /// @param[in] swing_start_time Start time of the reference tracking.
  /// @param[in] period_swing Period where the foot is swinging.
  /// @param[in] period_stance Period where the foot is stancing.
  /// @param[in] num_phases_in_period Number of phases in a period. Must be 
  /// positive. Default is 1.
  ///
  MPCPeriodicSwingFootRef(const int contact_index,
                          const double swing_height, const double swing_start_time, 
                          const double period_swing, const double period_stance,
                          const int num_phases_in_period=1);

  ///
  /// @brief Destructor. 
  ///
  ~MPCPeriodicSwingFootRef();

  DEFINE_DEFAULT_CLONE_TASK_SPACE_3D_REF(MPCPeriodicSwingFootRef)

  ///
  /// @brief Set period. 
  /// @param[in] swing_start_time Start time of the reference tracking.
  /// @param[in] period_swing Period where the leg is swinging.
  /// @param[in] period_stance Period where the leg is in stance.
  /// @param[in] num_phases_in_period Number of phases in a period. Must be 
  /// positive. Default is 1.
  ///
  void setPeriod(const double swing_start_time, const double period_swing, 
                 const double period_stance, const int num_phases_in_period=1);

  ///
  /// @brief Set the reference positions of CoM from the contact positions of 
  /// the contact sequence. Also, the CoM refs of the first and last contact 
  /// phases are defined by user.
  /// @param[in] contact_sequence Contact sequence.
  /// @param[in] foot_step_planner Foot step planner.
  ///
  void setSwingFootRef(const std::shared_ptr<ContactSequence>& contact_sequence,
                       const std::shared_ptr<ContactPlannerBase>& foot_step_planner);

  void updateRef(const GridInfo& grid_info, Eigen::VectorXd& x3d_ref) const override;

  bool isActive(const GridInfo& grid_info) const override;

private:
  int contact_index_;
  std::vector<Eigen::Vector3d> contact_position_;
  std::vector<bool> is_contact_active_;
  double swing_height_, swing_start_time_, period_swing_, period_stance_, period_;
  int num_phases_in_period_;
};

} // namespace robotoc

#endif // ROBOTOC_MPC_PERIODIC_SWING_FOOT_REF_HPP_ 