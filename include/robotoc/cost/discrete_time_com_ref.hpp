#ifndef ROBOTOC_DISCRETE_TIME_COM_REF_HPP_
#define ROBOTOC_DISCRETE_TIME_COM_REF_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/time_varying_com_cost.hpp"
#include "robotoc/hybrid/contact_sequence.hpp"


namespace robotoc {

///
/// @class DiscreteTimeCoMRef
/// @brief Discrete-time periodic reference of the center of mass. 
///
class DiscreteTimeCoMRef : public TimeVaryingCoMRefBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] com_ref0 Initial CoM position reference.
  /// @param[in] com_step Reference step of the CoM over a contact phase.
  /// @param[in] start_phase Start contact phase.
  /// @param[in] end_phase End contact phase.
  /// @param[in] active_phases Number of phases where the tracking is active.
  /// @param[in] inactive_phases Number of phases where the tracking is inactive.
  /// @param[in] is_first_move_half If true, the first reference CoM step is 
  /// half step. 
  ///
  DiscreteTimeCoMRef(
      const std::vector<Eigen::Vector3d>& com_position_to_foot_position);

  ///
  /// @brief Destructor. 
  ///
  ~DiscreteTimeCoMRef();

  ///
  /// @brief Set swing foot reference. 
  /// @param[in] contact_index Contact index of interest.
  /// @param[in] initial_contact_frame_position Contact frame position at the 
  /// initial time of the horizon. Used if the contact frame is in swinging 
  /// phase at the beginning of the horizon.
  ///
  void setCoMRef(
      const std::shared_ptr<ContactSequence>& contact_sequence, 
      const Eigen::Vector3d& initial_com_position=Eigen::Vector3d::Zero());

  void update_com_ref(const GridInfo& grid_info, 
                      Eigen::VectorXd& com_ref) const override;

  bool isActive(const GridInfo& grid_info) const override;

private:
  std::vector<Eigen::Vector3d> com_position_, com_position_to_foot_position_;
  Eigen::Vector3d initial_com_position_, com_avg_;
  double initial_rate_from_com_position_;
  std::vector<bool> has_active_contacts_, has_inactive_contacts_;

  int nextActiveContactPhase(const int contact_phase) const {
    for (int phase=contact_phase+1; phase<has_active_contacts_.size(); ++phase) {
      if (has_active_contacts_[phase]) {
        return phase;
      }
    }
    return has_active_contacts_.size();
  }

  static double planarDistance(const Eigen::Vector3d& a, 
                               const Eigen::Vector3d& b) {
    const double xdiff = a.coeff(0) - b.coeff(0);
    const double ydiff = a.coeff(1) - b.coeff(1);
    return std::sqrt(xdiff*xdiff + ydiff*ydiff);
  }

};

} // namespace robotoc

#endif // ROBOTOC_DISCRETE_TIME_COM_REF_HPP_