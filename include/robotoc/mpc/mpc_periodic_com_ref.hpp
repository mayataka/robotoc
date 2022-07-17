#ifndef ROBOTOC_MPC_PERIODIC_COM_REF_HPP_
#define ROBOTOC_MPC_PERIODIC_COM_REF_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/com_ref_base.hpp"
#include "robotoc/hybrid/contact_sequence.hpp"
#include "robotoc/mpc/contact_planner_base.hpp"


namespace robotoc {

///
/// @class MPCPeriodicCoMRef
/// @brief Periodic reference positions of the center of mass. 
///
class MPCPeriodicCoMRef final : public CoMRefBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] swing_start_time Start time of the reference tracking.
  /// @param[in] period_active Period where the tracking is active.
  /// @param[in] period_inactive Period where the tracking is inactive.
  /// @param[in] num_phases_in_period Number of phases in a period. Must be 
  /// positive. Default is 1.
  ///
  MPCPeriodicCoMRef(const double swing_start_time, const double period_active, 
                    const double period_inactive, 
                    const int num_phases_in_period=1);

  ///
  /// @brief Destructor. 
  ///
  ~MPCPeriodicCoMRef();

  DEFINE_DEFAULT_CLONE_COM_REF(MPCPeriodicCoMRef)

  ///
  /// @brief Set period. 
  /// @param[in] swing_start_time Start time of the reference tracking.
  /// @param[in] period_active Period where the tracking is active.
  /// @param[in] period_inactive Period where the tracking is inactive.
  /// @param[in] num_phases_in_period Number of phases in a period. Must be 
  /// positive. Default is 1.
  ///
  void setPeriod(const double swing_start_time, const double period_active, 
                 const double period_inactive, 
                 const int num_phases_in_period=1);

  ///
  /// @brief Set the reference positions of CoM from the contact positions of 
  /// the contact sequence. Also, the CoM refs of the first and last contact 
  /// phases are defined by user.
  /// @param[in] contact_sequence Contact sequence.
  /// @param[in] foot_step_planner Foot step planner.
  ///
  void setCoMRef(const std::shared_ptr<ContactSequence>& contact_sequence,
                 const std::shared_ptr<ContactPlannerBase>& foot_step_planner);

  void updateRef(const GridInfo& grid_info, 
                      Eigen::VectorXd& com_ref) const override;

  bool isActive(const GridInfo& grid_info) const override;

private:
  std::vector<Eigen::Vector3d> com_;
  std::vector<bool> has_inactive_contacts_;
  double swing_start_time_, period_active_, period_inactive_, period_;
  int num_phases_in_period_;
};

} // namespace robotoc

#endif // ROBOTOC_MPC_PERIODIC_COM_REF_HPP_