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
  /// @param[in] vcom_ref Reference veloity of the CoM.
  /// @param[in] t0 Start time of the reference tracking.
  /// @param[in] period_active Period where the tracking is active.
  /// @param[in] period_inactive Period where the tracking is inactive.
  /// @param[in] is_first_move_half If true, the first reference CoM movement is 
  /// half speed. 
  ///
  PeriodicCoMRef2(const Eigen::Vector3d com_ref0, const Eigen::Vector3d vcom_ref, 
                  const double t0, const double period_active, 
                  const double period_inactive, const bool is_first_move_half);

  ///
  /// @brief Destructor. 
  ///
  ~PeriodicCoMRef2();

  void update_com_ref(const double t, Eigen::VectorXd& com_ref) const override;

  bool isActive(const double t) const override;

private:
  Eigen::Vector3d com_ref0_, vcom_ref_;
  double step_length_, t0_, period_active_, period_inactive_, period_;
  bool is_first_move_half_;

};

} // namespace robotoc

#endif // ROBOTOC_PERIODIC_COM_REF2_HPP_ 