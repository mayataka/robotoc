#ifndef IDOCP_CONTACT_STATUS_HPP_
#define IDOCP_CONTACT_STATUS_HPP_

#include <vector>

#include "idocp/robot/robot.hpp"

namespace idocp {

class ContactStatus {
public:
  ContactStatus(const Robot& robot);

  ContactStatus();
 
  ~ContactStatus();

  ContactStatus(const ContactStatus&) = default;

  ContactStatus& operator=(const ContactStatus&) = default;

  ContactStatus(ContactStatus&&) noexcept = default;

  ContactStatus& operator=(ContactStatus&&) noexcept = default;

  bool isContactActive(const int contact_index) const;

  const std::vector<bool>& isContactActive() const;

  int dimf() const;

  bool hasActiveContacts() const;

  void setContactStatus(const std::vector<bool>& is_contact_active);

  void activateContact(const int contact_index);

  void deactivateContact(const int contact_index);

  void activateContacts(const std::vector<int>& contact_indices);

  void deactivateContacts(const std::vector<int>& contact_indices);

private:
  std::vector<bool> is_contact_active_;
  int dimf_, max_point_contacts_;
  bool has_active_contacts_;

  void set_has_active_contacts();

};

} // namespace idocp

#include "idocp/robot/contact_status.hxx"

#endif // IDOCP_CONTACT_STATUS_HPP_ 