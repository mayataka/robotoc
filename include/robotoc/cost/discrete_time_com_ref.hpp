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
/// @brief Discrete-time reference of the center of mass. 
///
class DiscreteTimeCoMRef : public TimeVaryingCoMRefBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] com_position_to_contact_position Relative contact positions 
  /// from the CoM position position.
  ///
  DiscreteTimeCoMRef(
      const std::vector<Eigen::Vector3d>& com_position_to_contact_position);

  ///
  /// @brief Destructor. 
  ///
  ~DiscreteTimeCoMRef();

  ///
  /// @brief Set the CoM reference positions from the contact positions of the
  /// contact sequence. The first and last contact phases must have active 
  /// contacts.
  /// @param[in] contact_sequence Contact sequence.
  ///
  void setCoMRef(const std::shared_ptr<ContactSequence>& contact_sequence);

  ///
  /// @brief Set the CoM reference positions from the contact positions of the
  /// contact sequence. Also, the CoM refs of the first and last contact phases 
  /// are defined by user.
  /// @param[in] contact_sequence Contact sequence.
  /// @param[in] first_com_ref CoM reference at the first contact phase.
  /// @param[in] last_com_ref CoM reference at the last contact phase.
  ///
  void setCoMRef(const std::shared_ptr<ContactSequence>& contact_sequence,
                 const Eigen::Vector3d& first_com_ref, 
                 const Eigen::Vector3d& last_com_ref, 
                 const double first_rate=0, const double last_rate=0);

  void update_com_ref(const GridInfo& grid_info, 
                      Eigen::VectorXd& com_ref) const override;

  bool isActive(const GridInfo& grid_info) const override;

private:
  std::vector<Eigen::Vector3d> com_position_, com_position_to_foot_position_;
  std::vector<bool> has_inactive_contacts_;
  double first_rate_, last_rate_;
  int num_contact_phases_;

};

} // namespace robotoc

#endif // ROBOTOC_DISCRETE_TIME_COM_REF_HPP_