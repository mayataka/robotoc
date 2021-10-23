#ifndef ROBOTOC_IMPULSE_STATUS_HXX_
#define ROBOTOC_IMPULSE_STATUS_HXX_

#include "robotoc/robot/impulse_status.hpp"

#include <cassert>
#include <random>


namespace robotoc {

inline ImpulseStatus::ImpulseStatus(const int max_point_contacts)
  : contact_status_(max_point_contacts) {
}


inline ImpulseStatus::ImpulseStatus() 
  : contact_status_() {
}
 

inline ImpulseStatus::~ImpulseStatus() {
}


inline bool ImpulseStatus::operator==(const ImpulseStatus& other) const {
  assert(other.maxPointContacts() == maxPointContacts());
  for (int i=0; i<maxPointContacts(); ++i) {
    if (other.isImpulseActive(i) != isImpulseActive(i)) {
      return false;
    }
    if (!other.contactPoints()[i].isApprox(contactPoints()[i])) {
      return false;
    }
  }
  return true;
}


inline bool ImpulseStatus::operator!=(const ImpulseStatus& other) const {
  return !(*this == other);
}


inline bool ImpulseStatus::isImpulseActive(const int contact_index) const {
  return contact_status_.isContactActive(contact_index);
}


inline const std::vector<bool>& ImpulseStatus::isImpulseActive() const {
  return contact_status_.isContactActive();
}


inline bool ImpulseStatus::hasActiveImpulse() const {
  return contact_status_.hasActiveContacts();
}


inline int ImpulseStatus::dimf() const {
  return contact_status_.dimf();
}


inline int ImpulseStatus::maxPointContacts() const {
  return contact_status_.maxPointContacts();
}


inline void ImpulseStatus::setActivity(
    const ContactStatus& pre_contact_status, 
    const ContactStatus& post_contact_status) {
  assert(pre_contact_status.maxPointContacts() == maxPointContacts());
  assert(post_contact_status.maxPointContacts() == maxPointContacts());
  for (int i=0; i<contact_status_.maxPointContacts(); ++i) {
    if (pre_contact_status.isContactActive(i)) {
      deactivateImpulse(i);
    }
    else {
      if (post_contact_status.isContactActive(i)) {
        activateImpulse(i);
      }
      else {
        deactivateImpulse(i);
      }
    }
  }
}


inline void ImpulseStatus::setActivity(
    const std::vector<bool>& is_impulse_active) {
  contact_status_.setActivity(is_impulse_active);
}


inline void ImpulseStatus::activateImpulse(const int impulse_index) {
  contact_status_.activateContact(impulse_index);
}


inline void ImpulseStatus::deactivateImpulse(const int impulse_index) {
  contact_status_.deactivateContact(impulse_index);
}


inline void ImpulseStatus::activateImpulses(
    const std::vector<int>& impulse_indices) {
  contact_status_.activateContacts(impulse_indices);
} 


inline void ImpulseStatus::activateImpulses() {
  contact_status_.activateContacts();
} 


inline void ImpulseStatus::deactivateImpulses(
    const std::vector<int>& impulse_indices) {
  contact_status_.deactivateContacts(impulse_indices);
}


inline void ImpulseStatus::deactivateImpulses() {
  contact_status_.deactivateContacts();
}


inline void ImpulseStatus::setContactPoint(
    const int contact_index, const Eigen::Vector3d& contact_point) {
  contact_status_.setContactPoint(contact_index, contact_point);
}


inline void ImpulseStatus::setContactPoints(
    const std::vector<Eigen::Vector3d>& contact_points) {
  contact_status_.setContactPoints(contact_points);
}


inline const Eigen::Vector3d& ImpulseStatus::contactPoint(
    const int impulse_index) const {
  return contact_status_.contactPoint(impulse_index);
}


inline const std::vector<Eigen::Vector3d>& 
ImpulseStatus::contactPoints() const {
  return contact_status_.contactPoints();
}


inline void ImpulseStatus::setRandom() {
  contact_status_.setRandom();
}

} // namespace robotoc

#endif // ROBOTOC_IMPULSE_STATUS_HXX_ 