#ifndef ROBOTOC_CONTACT_STATUS_HPP_
#define ROBOTOC_CONTACT_STATUS_HPP_

#include <vector>
#include <iostream>

#include "Eigen/Core"


namespace robotoc {

///
/// @enum ContactType 
/// @brief Types of contacts 
///
enum class ContactType {
  PointContact,
  SurfaceContact
};


///
/// @class ContactStatus
/// @brief Contact status of robot model.
///
class ContactStatus {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] max_point_contacts Maximum number of the point contacts. 
  /// @param[in] contact_id Identifier number of the contact. Can be used only 
  /// in user-defined cost and constraints. Default is 0.
  ///
  ContactStatus(const int max_point_contacts, const int contact_id=0);

  ///
  /// @brief Default constructor. 
  ///
  ContactStatus();

  ///
  /// @brief Destructor. 
  ///
  ~ContactStatus();

  ///
  /// @brief Default copy constructor. 
  ///
  ContactStatus(const ContactStatus&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  ContactStatus& operator=(const ContactStatus&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ContactStatus(ContactStatus&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ContactStatus& operator=(ContactStatus&&) noexcept = default;

  ///
  /// @brief Defines a comparison operator. 
  ///
  bool operator==(const ContactStatus& other) const;

  ///
  /// @brief Defines a comparison operator. 
  ///
  bool operator!=(const ContactStatus& other) const;

  ///
  /// @brief Returns true if a contact is active and false if not.
  /// @param[in] contact_index Index of a contact of interedted. 
  /// @return true if a contact is active and false if not. 
  ///
  bool isContactActive(const int contact_index) const;

  ///
  /// @brief Returns the activity of the contacts.
  /// @return Const reference to the activity of the contacts. 
  ///
  const std::vector<bool>& isContactActive() const;

  ///
  /// @brief Returns true if there are active contacts and false if not.
  /// @return true if there are active contacts and false if not. 
  ///
  bool hasActiveContacts() const;

  ///
  /// @brief Returns the dimension of the active contact forces.
  /// @return Dimension of the active contacts forces. 
  ///
  int dimf() const;

  ///
  /// @brief Returns the maximum number of the contacts.
  /// @return The maximum number of the contacts. 
  ///
  int maxPointContacts() const;

  ///
  /// @brief Sets the activity of the contacts.
  /// @param[in] is_contact_active Activity of the contacts. Size must be 
  /// ContactStatus::maxPointContacts();
  ///
  void setActivity(const std::vector<bool>& is_contact_active);

  ///
  /// @brief Activates a contact.
  /// @param[in] contact_index Index of the contact that is activated.
  ///
  void activateContact(const int contact_index);

  ///
  /// @brief Deactivates a contact.
  /// @param[in] contact_index Index of the contact that is deactivated.
  ///
  void deactivateContact(const int contact_index);

  ///
  /// @brief Activates contacts.
  /// @param[in] contact_indices Indices of the contacts that are activated.
  ///
  void activateContacts(const std::vector<int>& contact_indices);

  ///
  /// @brief Activates all contacts.
  ///
  void activateContacts();

  ///
  /// @brief Deactivates contacts.
  /// @param[in] contact_indices Indices of the contacts that are deactivated.
  ///
  void deactivateContacts(const std::vector<int>& contact_indices);

  ///
  /// @brief Deactivates all contacts.
  ///
  void deactivateContacts();

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
  /// ContactStatus::maxPointContacts().
  ///
  void setContactPoints(const std::vector<Eigen::Vector3d>& contact_points);

  ///
  /// @brief Gets contact point.
  /// @param[in] contact_index Index of the contact .
  /// @return const reference to the contact point. 
  ///
  const Eigen::Vector3d& contactPoint(const int contact_index) const;

  ///
  /// @brief Gets contact points.
  /// @return const reference to the vector of contact points. 
  ///
  const std::vector<Eigen::Vector3d>& contactPoints() const;

  ///
  /// @brief Sets a contact surface by its normal vector.
  /// @param[in] contact_index Index of the contact.
  /// @param[in] contact_surface_normal Normal vector of the contact surface.
  ///
  void setContactSurfaceNormal(const int contact_index, 
                               const Eigen::Vector3d& contact_surface_normal);

  ///
  /// @brief Sets contact surfaces by their normal vectors.
  /// @param[in] contact_surfaces_normals Normal vectors of the contact surfaces. 
  /// Size must be ContactStatus::maxPointContacts().
  ///
  void setContactSurfacesNormals(
      const std::vector<Eigen::Vector3d>& contact_surfaces_normals);

  ///
  /// @brief Gets normal vector of the contact surface.
  /// @param[in] contact_index Index of the contact .
  /// @return const reference to the normal vector of the contact surface. 
  ///
  const Eigen::Vector3d& contactSurfaceNormal(const int contact_index) const;

  ///
  /// @brief Gets normal vectors of the contact surfaces.
  /// @return const reference to the vector of the normal vectors of the contact 
  /// surfaces. 
  ///
  const std::vector<Eigen::Vector3d>& contactSurfacesNormals() const;

  ///
  /// @brief Sets contact id.
  /// @param[in] contact_id Contact id. 
  ///
  void setContactId(const int contact_id);

  ///
  /// @brief Gets contact id.
  /// @return Contact id. 
  ///
  int contactId() const;

  ///
  /// @brief Fills contact status randomly.
  ///
  void setRandom();

  ///
  /// @brief Displays the contact status onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const ContactStatus& contact_status);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::vector<bool> is_contact_active_;
  std::vector<Eigen::Vector3d> contact_points_;
  std::vector<Eigen::Vector3d> contact_surfaces_normals_;
  int dimf_, max_point_contacts_, contact_id_;
  bool has_active_contacts_;

  void set_has_active_contacts();

};

} // namespace robotoc

#include "robotoc/robot/contact_status.hxx"

#endif // ROBOTOC_CONTACT_STATUS_HPP_ 