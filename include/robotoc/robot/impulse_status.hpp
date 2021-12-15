#ifndef ROBOTOC_IMPULSE_STATUS_HPP_
#define ROBOTOC_IMPULSE_STATUS_HPP_

#include <vector>
#include <iostream>

#include "Eigen/Core"

#include "robotoc/robot/contact_status.hpp"


namespace robotoc {

  ///
  /// @class ImpulseStatus
  /// @brief Impulse status of robot model. Wrapper of ContactStatus to treat
  /// impulses.
  ///
class ImpulseStatus {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] max_point_contacts Maximum number of the point contacts. 
  /// @param[in] impulse_id Identifier number of the impulse. Can be used only 
  /// in user-defined cost and constraints. Default is 0.
  ///
  ImpulseStatus(const int max_point_contacts, const int impulse_id=0);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseStatus();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseStatus();

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseStatus(const ImpulseStatus&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  ImpulseStatus& operator=(const ImpulseStatus&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseStatus(ImpulseStatus&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseStatus& operator=(ImpulseStatus&&) noexcept = default;

  ///
  /// @brief Defines a comparison operator. 
  ///
  bool operator==(const ImpulseStatus& other) const;

  ///
  /// @brief Defines a comparison operator. 
  ///
  bool operator!=(const ImpulseStatus& other) const;

  ///
  /// @brief Returns true if a contact is active and false if not.
  /// @param[in] contact_index Index of a contact of interedted. 
  /// @return true if a contact is active and false if not. 
  ///
  bool isImpulseActive(const int contact_index) const;

  ///
  /// @brief Returns the activity of the impulses.
  /// @return Const reference to the activity of the impulses. 
  ///
  const std::vector<bool>& isImpulseActive() const;

  ///
  /// @brief Returns true if there are active impulses and false if not.
  /// @return true if there are active impulses and false if not. 
  ///
  bool hasActiveImpulse() const;

  ///
  /// @brief Returns the dimension of the active impulse forces.
  /// @return Dimension of the active impulse forces.  
  ///
  int dimi() const;

  ///
  /// @brief Returns the maximum number of the contacts.
  /// @return The maximum number of the contacts. 
  ///
  int maxPointContacts() const;

  ///
  /// @brief Sets the activity of the impulses from two sequential contact 
  /// statuses.
  /// @param[in] pre_contact_status Contact status before this impulse. 
  /// @param[in] post_contact_status Contact status after this impulse. 
  ///
  void setActivity(const ContactStatus& pre_contact_status, 
                   const ContactStatus& post_contact_status);

  ///
  /// @brief Sets the activity of the impulses.
  /// @param[in] is_impulse_active Activity of the impulses. Size must be 
  /// ContactStatus::maxPointContacts();
  ///
  void setActivity(const std::vector<bool>& is_impulse_active);

  ///
  /// @brief Activates a impulse.
  /// @param[in] contact_index Index of the contact that is activated.
  ///
  void activateImpulse(const int contact_index);

  ///
  /// @brief Deactivates a impulse.
  /// @param[in] contact_index Index of the contact that is deactivated.
  ///
  void deactivateImpulse(const int contact_index);

  ///
  /// @brief Activates impulses.
  /// @param[in] impulse_indices Indices of the impulses that are activated.
  ///
  void activateImpulses(const std::vector<int>& impulse_indices);

  ///
  /// @brief Activates all impulses.
  ///
  void activateImpulses();

  ///
  /// @brief Deactivates impulses.
  /// @param[in] impulse_indices Indices of the impulses that are deactivated.
  ///
  void deactivateImpulses(const std::vector<int>& impulse_indices);

  ///
  /// @brief Deactivates all impulses.
  ///
  void deactivateImpulses();

  ///
  /// @brief Sets a contact point.
  /// @param[in] contact_index Index of the contact.
  /// @param[in] contact_point Contact point.
  ///
  void setContactPoint(const int contact_index, 
                       const Eigen::Vector3d& contact_point);

  ///
  /// @brief Sets contact points.
  /// @param[in] contact_points Contact points. Size must be 
  /// ImpulseStatus::maxPointContacts().
  ///
  void setContactPoints(const std::vector<Eigen::Vector3d>& contact_points);

  ///
  /// @brief Gets contact point of impulse.
  /// @param[in] impulse_index Index of the impulse.
  /// @return const reference to the contact points. 
  ///
  const Eigen::Vector3d& contactPoint(const int impulse_index) const;

  ///
  /// @brief Gets contact points.
  /// @return const reference to the contact points. 
  ///
  const std::vector<Eigen::Vector3d>& contactPoints() const;

  ///
  /// @brief Sets the rotation matrix of a contact surface.
  /// @param[in] contact_index Index of the contact.
  /// @param[in] contact_surface_rotation Rotation matrix of the contact surface.
  ///
  void setContactSurfaceRotation(const int contact_index, 
                                 const Eigen::Matrix3d& contact_surface_rotation);

  ///
  /// @brief Sets the rotation matrices of contact surfaces.
  /// @param[in] contact_surfaces_rotations Rotation matrices of the contact 
  //// surfaces. Size must be ImpulseStatus::maxPointContacts().
  ///
  void setContactSurfacesRotations(
      const std::vector<Eigen::Matrix3d>& contact_surfaces_rotations);

  ///
  /// @brief Gets rotation matrix of a contact surface.
  /// @param[in] contact_index Index of the contact.
  /// @return const reference to the rotation matrix of the contact surface. 
  ///
  const Eigen::Matrix3d& contactSurfaceRotation(const int contact_index) const;

  ///
  /// @brief Gets rotation matrices of the contact surfaces.
  /// @return const reference to the rotation matrices of the contact surfaces. 
  ///
  const std::vector<Eigen::Matrix3d>& contactSurfacesRotations() const;

  ///
  /// @brief Sets impulse id.
  /// @param[in] impulse_id Impulse id. 
  ///
  void setImpulseId(const int impulse_id);

  ///
  /// @brief Gets impulse id.
  /// @return Impulse id. 
  ///
  int impulseId() const;

  ///
  /// @brief Fills impulse status randomly.
  ///
  void setRandom();

  ///
  /// @brief Displays the impulse status onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const ImpulseStatus& impulse_status);

private:
  ContactStatus contact_status_;

};

} // namespace robotoc

#include "robotoc/robot/impulse_status.hxx"

#endif // ROBOTOC_IMPULSE_STATUS_HPP_ 