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
  /// @brief Set the reference contact positions from the contact sequence. 
  /// The first and last contact phases must have active contacts.
  /// @param[in] contact_sequence Contact sequence.
  ///
  void setSwingFootRef(const std::shared_ptr<ContactSequence>& contact_sequence);

  ///
  /// @brief Set the reference contact positions from the contact sequence. 
  /// The first and last contact phases must have active contacts.
  /// Also, the position refs of the first and last contact phases are defined 
  /// by user.
  /// @param[in] contact_sequence Contact sequence.
  /// @param[in] first_contact_position Reference contact position at the first 
  /// contact phase.
  /// @param[in] last_contact_position Reference contact position at the last 
  /// contact phase.
  /// @param[in] first_rate Rate of the first contact position (0 <= rate <= 1).
  /// Default is 0.
  /// @param[in] last_rate Rate of the last contact position (0 <= rate <= 1).
  /// Default is 0.
  ///
  void setSwingFootRef(const std::shared_ptr<ContactSequence>& contact_sequence, 
                       const Eigen::Vector3d& first_contact_position, 
                       const Eigen::Vector3d& last_contact_position, 
                       const double first_rate=0, const double last_rate=0);

  void update_x3d_ref(const GridInfo& grid_info, 
                      Eigen::VectorXd& x3d_ref) const override;

  bool isActive(const GridInfo& grid_info) const override;

private:
  int contact_index_, num_contact_phases_;
  double step_height_, first_rate_, last_rate_;
  std::vector<Eigen::Vector3d> contact_position_;
  std::vector<bool> is_contact_active_;

};

} // namespace robotoc

#endif // ROBOTOC_DISCRETE_TIME_SWING_FOOT_REF_HPP_ 