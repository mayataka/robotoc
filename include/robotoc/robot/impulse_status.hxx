#ifndef ROBOTOC_IMPULSE_STATUS_HXX_
#define ROBOTOC_IMPULSE_STATUS_HXX_

#include "robotoc/robot/impulse_status.hpp"

#include <cassert>
#include <random>


namespace robotoc {

inline ImpulseStatus::ImpulseStatus(
    const std::vector<ContactType>& contact_types, const int impulse_mode_id)
  : contact_status_(contact_types, impulse_mode_id) {
}


inline ImpulseStatus::ImpulseStatus() 
  : contact_status_() {
}
 

inline ImpulseStatus::~ImpulseStatus() {
}


inline bool ImpulseStatus::operator==(const ImpulseStatus& other) const {
  assert(other.maxNumContacts() == maxNumContacts());
  for (int i=0; i<maxNumContacts(); ++i) {
    if (other.isImpulseActive(i) != isImpulseActive(i)) {
      return false;
    }
    if (!other.contactPlacement(i).isApprox(contactPlacement(i))) {
      return false;
    }
  }
  return true;
}


inline bool ImpulseStatus::operator!=(const ImpulseStatus& other) const {
  return !(*this == other);
}


inline ContactType ImpulseStatus::contactType(const int contact_index) const {
  return contact_status_.contactType(contact_index);
}


inline const std::vector<ContactType>& ImpulseStatus::contactTypes() const {
  return contact_status_.contactTypes();
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


inline int ImpulseStatus::dimi() const {
  return contact_status_.dimf();
}


inline int ImpulseStatus::maxNumContacts() const {
  return contact_status_.maxNumContacts();
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


inline void ImpulseStatus::deactivateImpulses(
    const std::vector<int>& impulse_indices) {
  contact_status_.deactivateContacts(impulse_indices);
}


inline void ImpulseStatus::setContactPlacement(
    const int contact_index, const Eigen::Vector3d& contact_position) {
  contact_status_.setContactPlacement(contact_index, contact_position);
}


inline void ImpulseStatus::setContactPlacement(
    const int contact_index, const Eigen::Vector3d& contact_position,
    const Eigen::Matrix3d& contact_rotation) {
  contact_status_.setContactPlacement(contact_index, contact_position, 
                                      contact_rotation);
}


inline void ImpulseStatus::setContactPlacement(const int contact_index, 
                                               const SE3& contact_placement) {
  contact_status_.setContactPlacement(contact_index, contact_placement);
}


inline void ImpulseStatus::setContactPlacements(
    const std::vector<Eigen::Vector3d>& contact_positions) {
  contact_status_.setContactPlacements(contact_positions);
}


inline void ImpulseStatus::setContactPlacements(
    const std::vector<Eigen::Vector3d>& contact_positions,
    const std::vector<Eigen::Matrix3d>& contact_rotations) {
  contact_status_.setContactPlacements(contact_positions, contact_rotations);
}


inline void ImpulseStatus::setContactPlacements(
    const aligned_vector<SE3>& contact_placements) {
  contact_status_.setContactPlacements(contact_placements);
}


inline const SE3& ImpulseStatus::contactPlacement(
    const int contact_index) const {
  return contact_status_.contactPlacement(contact_index);
}


inline const Eigen::Vector3d& ImpulseStatus::contactPosition(
    const int contact_index) const {
  return contact_status_.contactPosition(contact_index);
}


inline const Eigen::Matrix3d& ImpulseStatus::contactRotation(
    const int contact_index) const {
  return contact_status_.contactRotation(contact_index);
}


inline const aligned_vector<SE3>& ImpulseStatus::contactPlacements() const {
  return contact_status_.contactPlacements();
}


inline const std::vector<Eigen::Vector3d>& 
ImpulseStatus::contactPositions() const {
  return contact_status_.contactPositions();
}

 
inline const std::vector<Eigen::Matrix3d>& 
ImpulseStatus::contactRotations() const {
  return contact_status_.contactRotations();
}


inline void ImpulseStatus::setImpulseModeId(const int impulse_mode_id) {
  contact_status_.setContactModeId(impulse_mode_id);
}


inline int ImpulseStatus::impulseModeId() const {
  return contact_status_.contactModeId();
}


inline void ImpulseStatus::setRandom() {
  contact_status_.setRandom();
}

} // namespace robotoc

#endif // ROBOTOC_IMPULSE_STATUS_HXX_ 