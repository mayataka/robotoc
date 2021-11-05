#ifndef ROBOTOC_TROTTING_SWING_FOOT_REF_HPP_
#define ROBOTOC_TROTTING_SWING_FOOT_REF_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/swing_foot_cost.hpp"


namespace robotoc {

///
/// @class TrottingSwingFootRef
/// @brief Periodic reference of the foot position. 
///
class TrottingSwingFootRef : public SwingFootRefBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] x_ref_foot_contact_index Contact index of another foot that 
  /// gives the reference x coordinate.
  /// @param[in] y_ref_foot_contact_index Contact index of another foot that 
  /// gives the reference y coordinate.
  /// @param[in] step_length The step length of the gait.
  /// @param[in] step_height The step height of the gait.
  ///
  TrottingSwingFootRef(const int x_ref_foot_contact_index, 
                       const int y_ref_foot_contact_index,
                       const double step_length, const double step_height);

  ///
  /// @brief Destructor. 
  ///
  ~TrottingSwingFootRef();

  void update_q_3d_ref(const ContactStatus& contact_status, 
                       Eigen::VectorXd& q_3d_ref) const override;

private:
  int x_ref_foot_contact_index_, y_ref_foot_contact_index_;
  double step_length_, step_height_;

};

} // namespace robotoc

#endif // ROBOTOC_TROTTING_SWING_FOOT_REF_HPP_ 