#ifndef ROBOTOC_TROT_SWING_FOOT_REF_HPP_
#define ROBOTOC_TROT_SWING_FOOT_REF_HPP_

#include <string>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/swing_foot_cost.hpp"


namespace robotoc {

///
/// @class TrotSwingFootRef
/// @brief Periodic reference of the foot position. 
///
class TrotSwingFootRef : public SwingFootRefBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] contact_index Contact index of interest 
  /// @param[in] x_ref_foot_contact_index Contact index of another foot that 
  /// gives the reference x coordinate.
  /// @param[in] y_ref_foot_contact_index Contact index of another foot that 
  /// gives the reference y coordinate.
  /// @param[in] step_length The step length of the gait.
  /// @param[in] step_height The step height of the gait.
  ///
  TrotSwingFootRef(const int contact_index, 
                   const int x_ref_foot_contact_index, 
                   const int y_ref_foot_contact_index,
                   const double step_length, const double step_height);

  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] contact_frame_name Name of the contact frame of interest 
  /// @param[in] x_ref_foot_contact_index Contact index of another foot that 
  /// gives the reference x coordinate.
  /// @param[in] y_ref_foot_contact_index Contact index of another foot that 
  /// gives the reference y coordinate.
  /// @param[in] step_length The step length of the gait.
  /// @param[in] step_height The step height of the gait.
  ///
  TrotSwingFootRef(const Robot& robot, const std::string& contact_frame_name, 
                   const int x_ref_foot_contact_index, 
                   const int y_ref_foot_contact_index,
                   const double step_length, const double step_height);

  ///
  /// @brief Destructor. 
  ///
  ~TrotSwingFootRef();

  void update_x3d_ref(const ContactStatus& contact_status, 
                      Eigen::VectorXd& x3d_ref) const override;

private:
  int contact_index_, x_ref_foot_contact_index_, y_ref_foot_contact_index_;
  double step_length_, step_height_;

};

} // namespace robotoc

#endif // ROBOTOC_TROT_SWING_FOOT_REF_HPP_