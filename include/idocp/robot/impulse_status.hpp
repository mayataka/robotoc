#ifndef IDOCP_IMPULSE_STATUS_HPP_
#define IDOCP_IMPULSE_STATUS_HPP_

#include <vector>
#include "idocp/robot/contact_status.hpp"


namespace idocp {

  ///
  /// @class ImpulseStatus
  /// @brief Impulse status. A wrappper of ContactStatus for impulse condition.
  ///
class ImpulseStatus {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] max_point_contacts Maximum number of the point contacts. 
  ///
  ImpulseStatus(const int max_point_contacts);

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
  /// @brief Return true if a contact is active and false if not.
  /// @param[in] contact_index Index of a contact of interedted. 
  /// @return true if a contact is active and false if not. 
  ///
  bool isImpulseActive(const int contact_index) const;

  ///
  /// @brief Return contact status.
  /// @return Const reference to the contact status. 
  ///
  const std::vector<bool>& isImpulseActive() const;

  ///
  /// @brief Return true if there are active contacts and false if not.
  /// @return true if there are active contacts and false if not. 
  ///
  bool hasActiveImpulse() const;

  ///
  /// @brief Return the dimension of the active contacts.
  /// @return Dimension of the active contacts. 
  ///
  int dimp() const;

  ///
  /// @brief Return the number of the active contacts.
  /// @return The number of the active contacts. 
  ///
  int num_active_impulse() const;

  ///
  /// @brief Return the maximum number of the contacts.
  /// @return The maximum number of the contacts. 
  ///
  int max_point_contacts() const;

  ///
  /// @brief Set the contact status.
  /// @param[in] is_contact_active Contact status. Size must be 
  /// ImpulseStatus::max_point_contacts();
  ///
  void setImpulseStatus(const std::vector<bool>& is_impulse_active);

  ///
  /// @brief Activate a contact.
  /// @param[in] contact_index Index of the contact that is activated.
  ///
  void activateImpulse(const int contact_index);

  ///
  /// @brief Deactivate a contact.
  /// @param[in] contact_index Index of the contact that is deactivated.
  ///
  void deactivateImpulse(const int contact_index);

  ///
  /// @brief Activate contacts.
  /// @param[in] contact_indices Indices of the contacts that are activated.
  ///
  void activateImpulse(const std::vector<int>& contact_indices);

  ///
  /// @brief Deactivate contacts.
  /// @param[in] contact_indices Indices of the contacts that are deactivated.
  ///
  void deactivateImpulse(const std::vector<int>& contact_indices);

private:
  ContactStatus impulse_status_;

  void set_has_active_impulse();

};

} // namespace idocp

#include "idocp/robot/impulse_status.hxx"

#endif // IDOCP_IMPULSE_STATUS_HPP_ 