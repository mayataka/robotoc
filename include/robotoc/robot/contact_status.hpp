#ifndef ROBOTOC_CONTACT_STATUS_HPP_
#define ROBOTOC_CONTACT_STATUS_HPP_

#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>

#include "Eigen/Core"

#include "robotoc/robot/se3.hpp"
#include "robotoc/utils/aligned_vector.hpp"
#include "robotoc/utils/aligned_unordered_map.hpp"


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
  /// @param[in] contact_types Types of contacts. 
  /// @param[in] contact_frame_names Names of contact frames. 
  /// @param[in] default_friction_coefficient Default friction coefficitn. 
  /// Must be positive. Default is 0.7.
  /// @param[in] contact_mode_id Identifier number of the contact mode. Can be  
  /// used only in user-defined cost and constraints. Default is 0.
  ///
  ContactStatus(const std::vector<ContactType>& contact_types, 
                const std::vector<std::string>& contact_frame_names,
                const double default_friction_coefficient=0.7,
                const int contact_mode_id=0);

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
  /// @brief Returns the type of the contact.
  /// @param[in] contact_index Index of a contact of interedted. 
  /// @return Contact type. 
  ///
  ContactType contactType(const int contact_index) const;

  ///
  /// @brief Returns the types of the contacts.
  /// @return Contact types. 
  ///
  const std::vector<ContactType>& contactTypes() const;

  ///
  /// @brief Returns the name of the contact frame.
  /// @param[in] contact_index Index of a contact of interedted. 
  /// @return Name of the contact frame. 
  ///
  const std::string& contactFrameName(const int contact_index) const;

  ///
  /// @brief Returns the names of the contact frames.
  /// @return Name of the contact frames. 
  ///
  const std::vector<std::string>& contactFrameNames() const;

  ///
  /// @brief Returns true if a contact is active and false if not.
  /// @param[in] contact_index Index of a contact of interedted. 
  /// @return true if a contact is active and false if not. 
  ///
  bool isContactActive(const int contact_index) const;

  ///
  /// @brief Returns true if a contact is active and false if not.
  /// @param[in] contact_frame_name Name of the contact frame that is activated.
  /// @return true if a contact is active and false if not. 
  ///
  bool isContactActive(const std::string& contact_frame_name) const;

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
  int maxNumContacts() const;

  ///
  /// @brief Activates a contact.
  /// @param[in] contact_index Index of the contact that is activated.
  ///
  void activateContact(const int contact_index);

  ///
  /// @brief Activates a contact.
  /// @param[in] contact_frame_name Name of the contact frame that is activated.
  ///
  void activateContact(const std::string& contact_frame_name);

  ///
  /// @brief Deactivates a contact.
  /// @param[in] contact_index Index of the contact that is deactivated.
  ///
  void deactivateContact(const int contact_index);

  ///
  /// @brief Deactivates a contact.
  /// @param[in] contact_frame_name Name of the contact frame that is deactivated.
  ///
  void deactivateContact(const std::string& contact_frame_name);

  ///
  /// @brief Activates contacts.
  /// @param[in] contact_indices Indices of the contacts that are activated.
  ///
  void activateContacts(const std::vector<int>& contact_indices);

  ///
  /// @brief Activates contacts.
  /// @param[in] contact_frame_names Frame names of the contacts that are activated.
  ///
  void activateContacts(const std::vector<std::string>& contact_frame_names);

  ///
  /// @brief Deactivates contacts.
  /// @param[in] contact_indices Indices of the contacts that are deactivated.
  ///
  void deactivateContacts(const std::vector<int>& contact_indices);

  ///
  /// @brief Deactivates contacts.
  /// @param[in] contact_frame_names Frame names of the contacts that are deactivated.
  ///
  void deactivateContacts(const std::vector<std::string>& contact_frame_names);

  ///
  /// @brief Sets a contact placement, that is, the position and rotation of 
  /// the contact. The contact rotation is set to Eigen::Matrix3d::Identity(), 
  /// which represents the vertical direction to the ground. For the point 
  /// contacts, the rotation is only used in the friction cone constraints.
  /// For the surface contacts, the rotation represents the rotational contact
  /// constraints on the contact frame of the robot.
  /// @param[in] contact_index Index of the contact.
  /// @param[in] contact_position Contact position.
  ///
  void setContactPlacement(const int contact_index, 
                           const Eigen::Vector3d& contact_position);

  ///
  /// @brief Sets a contact placement, that is, the position and rotation of 
  /// the contact. The contact rotation is set to Eigen::Matrix3d::Identity(), 
  /// which represents the vertical direction to the ground. For the point 
  /// contacts, the rotation is only used in the friction cone constraints.
  /// For the surface contacts, the rotation represents the rotational contact
  /// constraints on the contact frame of the robot.
  /// @param[in] contact_frame_name Name of the contact frame.
  /// @param[in] contact_position Contact position.
  ///
  void setContactPlacement(const std::string& contact_frame_name, 
                           const Eigen::Vector3d& contact_position);

  ///
  /// @brief Sets a contact placement, that is, the position and rotation of 
  /// the contact. For the point contacts, the rotation is only used in the 
  /// friction cone constraints.
  /// For the surface contacts, the rotation represents the rotational contact
  /// constraints on the contact frame of the robot.
  /// @param[in] contact_index Index of the contact.
  /// @param[in] contact_position Contact position.
  /// @param[in] contact_rotation Contact rotation.
  ///
  void setContactPlacement(const int contact_index, 
                           const Eigen::Vector3d& contact_position, 
                           const Eigen::Matrix3d& contact_rotation);

  ///
  /// @brief Sets a contact placement, that is, the position and rotation of 
  /// the contact. For the point contacts, the rotation is only used in the 
  /// friction cone constraints.
  /// For the surface contacts, the rotation represents the rotational contact
  /// constraints on the contact frame of the robot.
  /// @param[in] contact_frame_name Name of the contact frame.
  /// @param[in] contact_position Contact position.
  /// @param[in] contact_rotation Contact rotation.
  ///
  void setContactPlacement(const std::string& contact_frame_name, 
                           const Eigen::Vector3d& contact_position, 
                           const Eigen::Matrix3d& contact_rotation);

  ///
  /// @brief Sets a contact placement, that is, the position and rotation of 
  /// the contact. For the point contacts, the rotation is only used in the 
  /// friction cone constraints.
  /// For the surface contacts, the rotation represents the rotational contact
  /// constraints on the contact frame of the robot.
  /// @param[in] contact_index Index of the contact.
  /// @param[in] contact_placement Contact placement.
  ///
  void setContactPlacement(const int contact_index, 
                           const SE3& contact_placement);

  ///
  /// @brief Sets a contact placement, that is, the position and rotation of 
  /// the contact. For the point contacts, the rotation is only used in the 
  /// friction cone constraints.
  /// For the surface contacts, the rotation represents the rotational contact
  /// constraints on the contact frame of the robot.
  /// @param[in] contact_frame_name Name of the contact frame.
  /// @param[in] contact_placement Contact placement.
  ///
  void setContactPlacement(const std::string& contact_frame_name, 
                           const SE3& contact_placement);

  ///
  /// @brief Sets contact placements. The rotation of each contact is set to
  /// Eigen::Matrix3d::Identity(), which represents the vertical direction
  /// to the ground.
  /// @param[in] contact_positions Contact positions. Size must be 
  /// ContactStatus::maxNumContacts().
  ///
  void setContactPlacements(
      const std::vector<Eigen::Vector3d>& contact_positions);

  ///
  /// @brief Sets contact placements. The rotation of each contact is set to
  /// Eigen::Matrix3d::Identity(), which represents the vertical direction
  /// to the ground.
  /// @param[in] contact_positions Contact positions. Size must be 
  /// ContactStatus::maxNumContacts().
  ///
  void setContactPlacements(
      const std::unordered_map<std::string, Eigen::Vector3d>& contact_positions);

  ///
  /// @brief Sets contact placements.
  /// @param[in] contact_positions Contact positions. Size must be 
  /// ContactStatus::maxNumContacts().
  /// @param[in] contact_rotations Contact rotations. Size must be 
  /// ContactStatus::maxNumContacts().
  ///
  void setContactPlacements(
      const std::vector<Eigen::Vector3d>& contact_positions,
      const std::vector<Eigen::Matrix3d>& contact_rotations);

  ///
  /// @brief Sets contact placements.
  /// @param[in] contact_positions Contact positions. Size must be 
  /// ContactStatus::maxNumContacts().
  /// @param[in] contact_rotations Contact rotations. Size must be 
  /// ContactStatus::maxNumContacts().
  ///
  void setContactPlacements(
      const std::unordered_map<std::string, Eigen::Vector3d>& contact_positions,
      const std::unordered_map<std::string, Eigen::Matrix3d>& contact_rotations);

  ///
  /// @brief Sets contact placements.
  /// @param[in] contact_placements Contact placements. Size must be 
  /// ContactStatus::maxNumContacts().
  ///
  void setContactPlacements(const aligned_vector<SE3>& contact_placements);
 
  ///
  /// @brief Sets contact placements.
  /// @param[in] contact_placements Contact placements. Size must be 
  /// ContactStatus::maxNumContacts().
  ///
  void setContactPlacements(
      const aligned_unordered_map<std::string, SE3>& contact_placements);

  ///
  /// @brief Gets the contact placement.
  /// @param[in] contact_index Index of the contact .
  /// @return const reference to the contact placement. 
  ///
  const SE3& contactPlacement(const int contact_index) const;

  ///
  /// @brief Gets the contact placement.
  /// @param[in] contact_frame_name Name of the contact frame.
  /// @return const reference to the contact placement. 
  ///
  const SE3& contactPlacement(const std::string& contact_frame_name) const;

  ///
  /// @brief Gets the contact position.
  /// @param[in] contact_index Index of the contact .
  /// @return const reference to the contact position. 
  ///
  const Eigen::Vector3d& contactPosition(const int contact_index) const;

  ///
  /// @brief Gets the contact position.
  /// @param[in] contact_frame_name Name of the contact frame.
  /// @return const reference to the contact position. 
  ///
  const Eigen::Vector3d& contactPosition(
      const std::string& contact_frame_name) const;

  ///
  /// @brief Gets the contact rotation.
  /// @param[in] contact_index Index of the contact .
  /// @return const reference to the contact rotation. 
  ///
  const Eigen::Matrix3d& contactRotation(const int contact_index) const;

  ///
  /// @brief Gets the contact rotation.
  /// @param[in] contact_frame_name Name of the contact frame.
  /// @return const reference to the contact rotation. 
  ///
  const Eigen::Matrix3d& contactRotation(
      const std::string& contact_frame_name) const;

  ///
  /// @brief Gets the contact placements.
  /// @return const reference to the contact placements. 
  ///
  const aligned_vector<SE3>& contactPlacements() const;

  ///
  /// @brief Gets the contact positions.
  /// @return const reference to the contact positions. 
  ///
  const std::vector<Eigen::Vector3d>& contactPositions() const;

  ///
  /// @brief Gets the contact rotations.
  /// @return const reference to the contact rotations. 
  ///
  const std::vector<Eigen::Matrix3d>& contactRotations() const;

  ///
  /// @brief Gets the friction coefficint.
  /// @param[in] contact_index Index of the contact.
  /// @param[in] friction_coefficient Friction coefficient. Must be positive.
  ///
  void setFrictionCoefficient(const int contact_index, 
                              const double friction_coefficient);

  ///
  /// @brief Gets the friction coefficint.
  /// @param[in] contact_frame_name Name of the contact frame.
  /// @param[in] friction_coefficient Friction coefficient. Must be positive.
  ///
  void setFrictionCoefficient(const std::string& contact_frame_name, 
                              const double friction_coefficient);

  ///
  /// @brief Sets the friction coefficints.
  /// @param[in] friction_coefficients Friction coefficients. 
  /// Size must be ContactStatus::maxNumContacts() and each element must be positive.
  ///
  void setFrictionCoefficients(const std::vector<double>& friction_coefficients);

  ///
  /// @brief Sets the friction coefficints.
  /// @param[in] friction_coefficients Friction coefficients. 
  /// Size must be ContactStatus::maxNumContacts() and each element must be positive.
  ///
  void setFrictionCoefficients(
      const std::unordered_map<std::string, double>& friction_coefficients);

  ///
  /// @brief Gets the friction coefficint. Default value is 0.7.
  /// @param[in] contact_index Index of the contact.
  /// @return Friction coefficient of the contact. 
  ///
  double frictionCoefficient(const int contact_index) const;

  ///
  /// @brief Gets the friction coefficint. Default value is 0.7.
  /// @param[in] contact_frame_name Name of the contact frame.
  /// @return Friction coefficient of the contact. 
  ///
  double frictionCoefficient(const std::string& contact_frame_name) const;

  ///
  /// @brief Gets the friction coefficints. Default value is 0.7.
  /// @return Friction coefficients of the contacts. 
  ///
  const std::vector<double>& frictionCoefficients() const;

  ///
  /// @brief Finds the contact index correspoinding to the input contact frame name.
  /// @param[in] contact_frame_name Name of the contact frame.
  /// @return Contact index. 
  ///
  int findContactIndex(const std::string& contact_frame_name) const;

  ///
  /// @brief Sets contact mode id.
  /// @param[in] contact_mode_id Contact mode id. 
  /// @note Default contact mode id is 0.
  ///
  void setContactModeId(const int contact_mode_id);

  ///
  /// @brief Gets contact mode id.
  /// @return Contact mode id. 
  /// @note Default contact mode id is 0.
  ///
  int contactModeId() const;

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
  std::vector<ContactType> contact_types_;
  std::vector<std::string> contact_frame_names_;
  std::vector<bool> is_contact_active_;
  aligned_vector<SE3> contact_placements_;
  std::vector<Eigen::Vector3d> contact_positions_;
  std::vector<Eigen::Matrix3d> contact_rotations_;
  std::vector<double> friction_coefficients_;
  int dimf_, max_contacts_, max_num_contacts_, contact_mode_id_;
  bool has_active_contacts_;

  void setHasActiveContacts();

};

} // namespace robotoc

#include "robotoc/robot/contact_status.hxx"

#endif // ROBOTOC_CONTACT_STATUS_HPP_ 