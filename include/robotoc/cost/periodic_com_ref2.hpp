#ifndef ROBOTOC_PERIODIC_COM_REF2_HPP_
#define ROBOTOC_PERIODIC_COM_REF2_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/time_varying_com_cost.hpp"


namespace robotoc {

///
/// @class PeriodicCoMRef2
/// @brief Periodic reference of the center of mass. 
///
class PeriodicCoMRef2 : public TimeVaryingCoMRefBase {
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
  PeriodicCoMRef2(const Eigen::Vector3d com_ref0, 
                  const Eigen::Vector3d com_step, 
                  const int start_phase, const int end_phase, 
                  const int active_phases, const int inactive_phases, 
                  const bool is_first_move_half);

  ///
  /// @brief Destructor. 
  ///
  ~PeriodicCoMRef2();

  void update_com_ref(const GridInfo& grid_info, 
                      Eigen::VectorXd& com_ref) const override;

  bool isActive(const GridInfo& grid_info) const override;

private:
  Eigen::Vector3d com_ref0_, com_step_;
  int start_phase_, end_phase_, active_phases_, inactive_phases_; 
  bool is_first_move_half_;

};

} // namespace robotoc

#endif // ROBOTOC_PERIODIC_COM_REF2_HPP_ 