#ifndef ROBOTOC_PERIODIC_FOOT_TRACK_REF2_HPP_
#define ROBOTOC_PERIODIC_FOOT_TRACK_REF2_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/time_varying_task_space_3d_cost.hpp"


namespace robotoc {

///
/// @class PeriodicFootTrackRef2
/// @brief Periodic reference of the foot position. 
///
class PeriodicFootTrackRef2 : public TimeVaryingTaskSpace3DRefBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] x3d0 Initial foot position reference.
  /// @param[in] step_length The step length of the gait.
  /// @param[in] step_height The step height of the gait.
  /// @param[in] start_phase Start contact phase.
  /// @param[in] end_phase End contact phase.
  /// @param[in] active_phases Number of phases where the tracking is active.
  /// @param[in] inactive_phases Number of phases where the tracking is inactive.
  /// @param[in] is_first_step_half If true, the length ofh te first reference 
  /// foot step is half. 
  ///
  PeriodicFootTrackRef2(const Eigen::Vector3d x3d0, 
                        const double step_length, const double step_height, 
                        const int start_phase, const int end_phase, 
                        const int active_phases, const int inactive_phases, 
                        const bool is_first_step_half);

  ///
  /// @brief Destructor. 
  ///
  ~PeriodicFootTrackRef2();

  void update_x3d_ref(const GridInfo& grid_info, 
                      Eigen::VectorXd& x3d_ref) const override;

  bool isActive(const GridInfo& grid_info) const override;

private:
  Eigen::Vector3d x3d0_;
  double step_length_, step_height_;
  int start_phase_, end_phase_, active_phases_, inactive_phases_; 
  bool is_first_step_half_;

};

} // namespace robotoc

#endif // ROBOTOC_PERIODIC_FOOT_TRACK_REF2_HPP_ 