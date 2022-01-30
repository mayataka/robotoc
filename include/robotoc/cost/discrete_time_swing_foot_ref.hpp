#ifndef ROBOTOC_DISCRETE_TIME_SWING_FOOT_REF_HPP_
#define ROBOTOC_DISCRETE_TIME_SWING_FOOT_REF_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/time_varying_task_space_3d_cost.hpp"
#include "robotoc/hybrid/contact_sequence.hpp"


namespace robotoc {

///
/// @class DiscreteTimeSwingFootRef
/// @brief Discrete-time reference of the positions of swinging contact frames. 
///
class DiscreteTimeSwingFootRef : public TimeVaryingTaskSpace3DRefBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] contact_index Contact index of interest.
  /// @param[in] step_height The step height of the gait.
  ///
  DiscreteTimeSwingFootRef(const int contact_index, const double step_height);

  ///
  /// @brief Destructor. 
  ///
  ~DiscreteTimeSwingFootRef();

  ///
  /// @brief Set swing foot reference. 
  /// @param[in] contact_index Contact index of interest.
  /// @param[in] initial_contact_frame_position Contact frame position at the 
  /// initial time of the horizon. Used if the contact frame is in swinging 
  /// phase at the beginning of the horizon.
  /// @param[in] contact_position_before_initial_time Contact frame position 
  /// before the initial time of the horizon. Used if the initial contact phase 
  /// is inactive.
  ///
  void setSwingFootRef(
      const std::shared_ptr<ContactSequence>& contact_sequence, 
      const Eigen::Vector3d& initial_contact_frame_position=Eigen::Vector3d::Zero(), 
      const Eigen::Vector3d& contact_position_before_initial_time=Eigen::Vector3d::Zero());

  void update_x3d_ref(const GridInfo& grid_info, 
                      Eigen::VectorXd& x3d_ref) const override;

  bool isActive(const GridInfo& grid_info) const override;

private:
  std::vector<Eigen::Vector3d> contact_position_;
  Eigen::Vector3d initial_contact_frame_position_;
  double step_height_, initial_rate_from_step_length_;
  int contact_index_;
  std::vector<bool> is_contact_active_;

  int nextActiveContactPhase(const int contact_phase) const {
    for (int phase=contact_phase+1; phase<is_contact_active_.size(); ++phase) {
      if (is_contact_active_[phase]) {
        return phase;
      }
    }
    return is_contact_active_.size();
  }

  static double planarDistance(const Eigen::Vector3d& a, 
                               const Eigen::Vector3d& b) {
    const double xdiff = a.coeff(0) - b.coeff(0);
    const double ydiff = a.coeff(1) - b.coeff(1);
    return std::sqrt(xdiff*xdiff + ydiff*ydiff);
  }

};

} // namespace robotoc

#endif // ROBOTOC_DISCRETE_TIME_SWING_FOOT_REF_HPP_ 