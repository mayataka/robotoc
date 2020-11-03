#ifndef IDOCP_IMPULSE_STATUS_HXX_
#define IDOCP_IMPULSE_STATUS_HXX_

#include "idocp/robot/impulse_status.hpp"

namespace idocp {

inline ImpulseStatus::ImpulseStatus(const int max_point_contacts)
  : impulse_status_(max_point_contacts) {
}


inline ImpulseStatus::ImpulseStatus() 
  : impulse_status_() {
}
 

inline ImpulseStatus::~ImpulseStatus() {
}


inline bool ImpulseStatus::isImpulseActive(const int contact_index) const {
  return impulse_status_.isContactActive(contact_index);
}


inline const std::vector<bool>& ImpulseStatus::isImpulseActive() const {
  return impulse_status_.isContactActive();
}


inline bool ImpulseStatus::hasActiveImpulse() const {
  return impulse_status_.hasActiveContacts();
}


inline int ImpulseStatus::dimp() const {
  return impulse_status_.dimf();
}


inline int ImpulseStatus::num_active_impulse() const {
  return impulse_status_.num_active_contacts();
}


inline int ImpulseStatus::max_point_contacts() const {
  return impulse_status_.max_point_contacts();
}


inline void ImpulseStatus::setImpulseStatus(
    const std::vector<bool>& is_impulse_active) {
  impulse_status_.setContactStatus(is_impulse_active);
}


inline void ImpulseStatus::activateImpulse(const int contact_index) {
  impulse_status_.activateContact(contact_index);
}


inline void ImpulseStatus::deactivateImpulse(const int contact_index) {
  impulse_status_.deactivateContact(contact_index);
}


inline void ImpulseStatus::activateImpulse(
    const std::vector<int>& contact_indices) {
  impulse_status_.activateContacts(contact_indices);
} 


inline void ImpulseStatus::deactivateImpulse(
    const std::vector<int>& contact_indices) {
  impulse_status_.deactivateContacts(contact_indices);
}

} // namespace idocp

#endif // IDOCP_IMPULSE_STATUS_HXX_ 