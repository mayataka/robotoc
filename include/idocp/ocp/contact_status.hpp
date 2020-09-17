#ifndef IDOCP_CONTACT_STATUS_HPP_
#define IDOCP_CONTACT_STATUS_HPP_

#include <vector>

#include "idocp/robot/robot.hpp"


namespace idocp {
class ContactStatus {
private:

  ContactStatus(const Robot& robot)
    : contact_status_(robot.max_point_contacts(), false) {
  }

  ~ContactStatus() {
  }

  const std::vector<bool>& ContactStatus() const {
    return contact_status_;
  }


  int num_active_contacts() const {
  }
  

  bool isContactActive(const int contact_index) const {
    assert(contact_index >= 0);
    assert(contact_index < contact_status_.size());
    return contact_status_[contact_index];
  }


  void setContactStatu(const std::vector<bool>& contact_status) {
    assert(contact_status.size() == contact_status_);
    contact_status_ = contact_status;
  }


  void activate(const int contact_index) {
    assert(contact_index >= 0);
    assert(contact_index < contact_status_.size());
    contact_status_[contact_index] = true;
  }


  void deactivate(const int contact_index) {
    assert(contact_index >= 0);
    assert(contact_index < contact_status_.size());
    contact_status_[contact_index] = false;
  }


public:
  std::vector<bool> contact_status_;
  std::vector<int> stack_contact_force_dim_;

};

} // namespace idocp 

#endif // IDOCP_CONTACT_STATUS_HPP_ 