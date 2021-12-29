#ifndef ROBOTOC_CONTACT_STATUS_HXX_
#define ROBOTOC_CONTACT_STATUS_HXX_

#include "robotoc/robot/contact_status.hpp"

#include <cassert>
#include <random>
#include <chrono>


namespace robotoc {

inline ContactStatus::ContactStatus(
    const std::vector<ContactType>& contact_types, const int contact_mode_id)
  : contact_types_(contact_types),
    is_contact_active_(contact_types.size(), false),
    contact_placements_(contact_types.size(), SE3::Identity()),
    contact_positions_(contact_types.size(), Eigen::Vector3d::Zero()),
    contact_rotations_(contact_types.size(), Eigen::Matrix3d::Identity()),
    dimf_(0),
    max_contacts_(contact_types.size()),
    max_num_contacts_(contact_types.size()),
    contact_mode_id_(contact_mode_id),
    has_active_contacts_(false) {
}


inline ContactStatus::ContactStatus() 
  : contact_types_(),
    is_contact_active_(),
    contact_placements_(),
    contact_positions_(),
    contact_rotations_(),
    dimf_(0),
    max_contacts_(0),
    max_num_contacts_(0),
    contact_mode_id_(0),
    has_active_contacts_(false) {
}


inline ContactStatus::~ContactStatus() {
}


inline bool ContactStatus::operator==(const ContactStatus& other) const {
  assert(other.maxNumContacts() == max_num_contacts_);
  for (int i=0; i<max_num_contacts_; ++i) {
    if (other.isContactActive(i) != isContactActive(i)) {
      return false;
    }
    if (!other.contactPlacement(i).isApprox(contactPlacement(i))) {
      return false;
    }
  }
  return true;
}


inline bool ContactStatus::operator!=(const ContactStatus& other) const {
  return !(*this == other);
}


inline ContactType ContactStatus::contactType(const int contact_index) const {
  assert(contact_index >= 0);
  assert(contact_index < max_num_contacts_);
  return contact_types_[contact_index];
}


inline const std::vector<ContactType>& ContactStatus::contactTypes() const {
  return contact_types_;
}


inline bool ContactStatus::isContactActive(const int contact_index) const {
  assert(!is_contact_active_.empty());
  assert(contact_index >= 0);
  assert(contact_index < is_contact_active_.size());
  return is_contact_active_[contact_index];
}


inline const std::vector<bool>& ContactStatus::isContactActive() const {
  return is_contact_active_;
}


inline bool ContactStatus::hasActiveContacts() const {
  return has_active_contacts_;
}


inline int ContactStatus::dimf() const {
  return dimf_;
}


inline int ContactStatus::maxNumContacts() const {
  return max_num_contacts_;
}


inline void ContactStatus::activateContact(const int contact_index) {
  assert(contact_index >= 0);
  assert(contact_index < max_num_contacts_);
  if (!is_contact_active_[contact_index]) {
    is_contact_active_[contact_index] = true;
    switch (contact_types_[contact_index]) {
      case ContactType::PointContact:
        dimf_ += 3;
        break;
      case ContactType::SurfaceContact:
        dimf_ += 6;
      default:
        break;
    }
  }
  set_has_active_contacts();
}


inline void ContactStatus::deactivateContact(const int contact_index) {
  assert(contact_index >= 0);
  assert(contact_index < max_num_contacts_);
  if (is_contact_active_[contact_index]) {
    is_contact_active_[contact_index] = false;
    switch (contact_types_[contact_index]) {
      case ContactType::PointContact:
        dimf_ -= 3;
        break;
      case ContactType::SurfaceContact:
        dimf_ -= 6;
      default:
        break;
    }
  }
  set_has_active_contacts();
}


inline void ContactStatus::activateContacts(
    const std::vector<int>& contact_indices) {
  assert(contact_indices.size() <= max_num_contacts_);
  for (const int e : contact_indices) {
    activateContact(e);
  }
}
 

inline void ContactStatus::deactivateContacts(
    const std::vector<int>& contact_indices) {
  assert(contact_indices.size() <= max_num_contacts_);
  for (const int e : contact_indices) {
    deactivateContact(e);
  }
}


inline void ContactStatus::setContactPlacement(
    const int contact_index, const Eigen::Vector3d& contact_position) {
  setContactPlacement(contact_index, contact_position, 
                      Eigen::Matrix3d::Identity());
}


inline void ContactStatus::setContactPlacement(
    const int contact_index, const Eigen::Vector3d& contact_position,
    const Eigen::Matrix3d& contact_rotation) {
  assert(contact_index >= 0);
  assert(contact_index < max_num_contacts_);
  contact_positions_[contact_index] = contact_position;
  contact_rotations_[contact_index] = contact_rotation;
  contact_placements_[contact_index] = SE3(contact_rotation, contact_position);
}


inline void ContactStatus::setContactPlacements(
    const std::vector<Eigen::Vector3d>& contact_positions) {
  assert(contact_positions.size() == max_num_contacts_);
  for (int i=0; i<max_num_contacts_; ++i) {
    setContactPlacement(i, contact_positions[i], Eigen::Matrix3d::Identity());
  }
}


inline void ContactStatus::setContactPlacements(
    const std::vector<Eigen::Vector3d>& contact_positions,
    const std::vector<Eigen::Matrix3d>& contact_rotations) {
  assert(contact_positions.size() == max_num_contacts_);
  assert(contact_rotations.size() == max_num_contacts_);
  for (int i=0; i<max_num_contacts_; ++i) {
    setContactPlacement(i, contact_positions[i], contact_rotations[i]);
  }
}


inline const SE3& ContactStatus::contactPlacement(
    const int contact_index) const {
  return contact_placements_[contact_index];
}


inline const Eigen::Vector3d& ContactStatus::contactPosition(
    const int contact_index) const {
  return contact_positions_[contact_index];
}


inline const Eigen::Matrix3d& ContactStatus::contactRotation(
    const int contact_index) const {
  return contact_rotations_[contact_index];
}


inline const aligned_vector<SE3>& ContactStatus::contactPlacements() const {
  return contact_placements_;
}


inline const std::vector<Eigen::Vector3d>& 
ContactStatus::contactPositions() const {
  return contact_positions_;
}


inline const std::vector<Eigen::Matrix3d>& 
ContactStatus::contactRotations() const {
  return contact_rotations_;
}


inline void ContactStatus::setContactModeId(const int contact_mode_id) {
  contact_mode_id_ = contact_mode_id;
}


inline int ContactStatus::contactModeId() const {
  return contact_mode_id_;
}


inline void ContactStatus::setRandom() {
  std::minstd_rand0 rand(
      std::chrono::system_clock::now().time_since_epoch().count());
  for (int i=0; i<max_num_contacts_; ++i) {
    if (rand()%2 == 0) {
      activateContact(i);
    }
    else {
      deactivateContact(i);
    }
  }
}


inline void ContactStatus::set_has_active_contacts() {
  has_active_contacts_ = (dimf_ > 0);
}

} // namespace robotoc

#endif // ROBOTOC_CONTACT_STATUS_HXX_ 