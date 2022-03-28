#ifndef ROBOTOC_MPC_PERIODIC_CONFIGURATION_REF_HPP_
#define ROBOTOC_MPC_PERIODIC_CONFIGURATION_REF_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/time_varying_configuration_space_cost.hpp"
#include "robotoc/hybrid/contact_sequence.hpp"
#include "robotoc/mpc/foot_step_planner_base.hpp"
#include "robotoc/utils/aligned_vector.hpp"


namespace robotoc {

///
/// @class MPCPeriodicConfigurationRef
/// @brief Periodic reference configuration for MPC of legged robots. 
///
class MPCPeriodicConfigurationRef final : public TimeVaryingConfigurationRefBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] q Reference configuration.
  /// @param[in] swing_start_time Start time of the reference tracking.
  /// @param[in] period_active Period where the tracking is active.
  /// @param[in] period_inactive Period where the tracking is inactive.
  ///
  MPCPeriodicConfigurationRef(const Eigen::VectorXd& q,
                              const double swing_start_time, 
                              const double period_active, 
                              const double period_inactive);

  ///
  /// @brief Destructor. 
  ///
  ~MPCPeriodicConfigurationRef();

  ///
  /// @brief Set period. 
  /// @param[in] swing_start_time Start time of the reference tracking.
  /// @param[in] period_active Period where the tracking is active.
  /// @param[in] period_inactive Period where the tracking is inactive.
  ///
  void setPeriod(const double swing_start_time, const double period_active, 
                 const double period_inactive);

  ///
  /// @brief Set the reference positions of CoM from the contact positions of 
  /// the contact sequence. Also, the CoM refs of the first and last contact 
  /// phases are defined by user.
  /// @param[in] contact_sequence Contact sequence.
  /// @param[in] foot_step_planner Foot step planner.
  ///
  void setConfigurationRef(const std::shared_ptr<ContactSequence>& contact_sequence,
                           const std::shared_ptr<FootStepPlannerBase>& foot_step_planner);

  void update_q_ref(const Robot& robot, const GridInfo& grid_info,
                    Eigen::VectorXd& q_ref) const override;

  bool isActive(const GridInfo& grid_info) const override;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::VectorXd q_;
  aligned_vector<Eigen::Quaterniond> quat_;
  std::vector<bool> has_inactive_contacts_;
  double swing_start_time_, period_active_, period_inactive_, period_;

};

} // namespace robotoc

#endif // RROBOTOC_MPC_PERIODIC_CONFIGURATION_REF_HPP_