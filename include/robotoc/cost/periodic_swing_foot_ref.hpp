#ifndef ROBOTOC_PERIODIC_SWING_FOOT_REF_HPP_
#define ROBOTOC_PERIODIC_SWING_FOOT_REF_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/task_space_3d_cost.hpp"


namespace robotoc {

///
/// @class PeriodicSwingFootRef
/// @brief Periodic reference of the swing foot position. 
///
class PeriodicSwingFootRef : public TaskSpace3DRefBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] x3d0 Initial foot position reference.
  /// @param[in] step_length The step length of the gait.
  /// @param[in] step_height The step height of the gait.
  /// @param[in] t0 Start time of the reference tracking.
  /// @param[in] period_swing Period where the foot is swinging.
  /// @param[in] period_stance Period where the foot is stancing.
  /// @param[in] is_first_step_half If true, the length ofh te first reference 
  /// foot step is half. 
  ///
  PeriodicSwingFootRef(const Eigen::Vector3d& x3d0, 
                       const Eigen::Vector3d& step_length, 
                       const double step_height, const double t0, 
                       const double period_swing, const double period_stance, 
                       const bool is_first_step_half);

  ///
  /// @brief Destructor. 
  ///
  ~PeriodicSwingFootRef();

  DEFINE_DEFAULT_CLONE_TASK_SPACE_3D_REF(PeriodicSwingFootRef)

  ///
  /// @brief Sets parameters. 
  /// @param[in] x3d0 Initial foot position reference.
  /// @param[in] step_length The step length of the gait.
  /// @param[in] step_height The step height of the gait.
  /// @param[in] t0 Start time of the reference tracking.
  /// @param[in] period_swing Period where the foot is swinging.
  /// @param[in] period_stance Period where the foot is stancing.
  /// @param[in] is_first_step_half If true, the length ofh te first reference 
  /// foot step is half. Default is false.
  ///
  void setFootTrackRef(const Eigen::Vector3d& x3d0, 
                       const Eigen::Vector3d& step_length, 
                       const double step_height, const double t0, 
                       const double period_swing, const double period_stance, 
                       const bool is_first_step_half=false);

  void updateRef(const GridInfo& grid_info, Eigen::VectorXd& x3d_ref) const override;

  bool isActive(const GridInfo& grid_info) const override;

private:
  Eigen::Vector3d x3d0_, step_length_;
  double step_height_, t0_, period_swing_, period_stance_, period_;
  bool is_first_step_half_;

};

} // namespace robotoc

#endif // ROBOTOC_PERIODIC_SWING_FOOT_REF_HPP_