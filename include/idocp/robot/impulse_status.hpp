#ifndef IDOCP_IMPULSE_STATUS_HPP_
#define IDOCP_IMPULSE_STATUS_HPP_

#include <vector>
#include "idocp/robot/contact_status.hpp"


namespace idocp {

  ///
  /// @class ImpulseStatus
  /// @brief Impulse status. A wrappper of ContactStatus for impulse problem.
  ///
class ImpulseStatus {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] max_point_contacts Maximum number of the point contacts. 
  ///
  ImpulseStatus(const int max_point_contacts);

  ///
  /// @brief Constructor. 
  /// @param[in] is_impulse_active Vector containing bool representing that each 
  /// point impulse is active or not. 
  ///
  ImpulseStatus(const std::vector<bool>& is_impulse_active);

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
  /// @brief Define comparison operator. 
  ///
  bool operator==(const ImpulseStatus& other) const;

  ///
  /// @brief Define comparison operator. 
  ///
  bool operator!=(const ImpulseStatus& other) const;

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
  /// @brief Return true if there are active impulse and false if not.
  /// @return true if there are active impulse and false if not. 
  ///
  bool hasActiveImpulse() const;

  ///
  /// @brief Return the dimension of the active impulse.
  /// @return Dimension of the active impulse. 
  ///
  int dimp() const;

  ///
  /// @brief Return the number of the active impulse.
  /// @return The number of the active impulse. 
  ///
  int num_active_impulse() const;

  ///
  /// @brief Return the maximum number of the contacts.
  /// @return The maximum number of the contacts. 
  ///
  int max_point_contacts() const;

  ///
  /// @brief Set the contact status from two sequential contact status.
  /// @param[in] current_contact_status Current contact status. 
  /// @param[in] next_contact_status Next contact status. 
  ///
  void setImpulseStatus(const ContactStatus& current_contact_status, 
                        const ContactStatus& next_contact_status);

  ///
  /// @brief Set the impulse status.
  /// @param[in] is_impulse_active Impulse status. Size must be 
  /// ImpulseStatus::max_point_contacts();
  ///
  void setImpulseStatus(const std::vector<bool>& is_impulse_active);

  ///
  /// @brief Activate a impulse.
  /// @param[in] contact_index Index of the contact that is activated.
  ///
  void activateImpulse(const int contact_index);

  ///
  /// @brief Deactivate a impulse.
  /// @param[in] contact_index Index of the contact that is deactivated.
  ///
  void deactivateImpulse(const int contact_index);

  ///
  /// @brief Activate impulse.
  /// @param[in] contact_indices Indices of the contacts that are activated.
  ///
  void activateImpulse(const std::vector<int>& contact_indices);

  ///
  /// @brief Deactivate impulse.
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