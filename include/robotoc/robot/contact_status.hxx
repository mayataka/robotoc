#ifndef ROBOTOC_CONTACT_STATUS_HXX_
#define ROBOTOC_CONTACT_STATUS_HXX_

#include "robotoc/robot/contact_status.hpp"

#include <stdexcept>
#include <cassert>
#include <random>
#include <chrono>


namespace robotoc {

inline ContactStatus::ContactStatus(
    const std::vector<ContactType>& contact_types, 
    const std::vector<std::string>& contact_frame_names,
    const double default_friction_coefficients)
  : contact_types_(contact_types),
    contact_frame_names_(contact_frame_names),
    is_contact_active_(contact_types.size(), false),
    contact_placements_(contact_types.size(), SE3::Identity()),
    contact_positions_(contact_types.size(), Eigen::Vector3d::Zero()),
    contact_rotations_(contact_types.size(), Eigen::Matrix3d::Identity()),
    friction_coefficients_(contact_types.size(), 0.7),
    dimf_(0),
    max_contacts_(contact_types.size()),
    max_num_contacts_(contact_types.size()),
    has_active_contacts_(false) {
  if (contact_types.size() != contact_frame_names.size()) {
    throw std::invalid_argument(
        "[ContactStatus] invalid argument: contact_frame_names.size() must be the same as contact_types.size()!");
  }
  if (default_friction_coefficients <= 0.0) {
    throw std::invalid_argument("[ContactStatus] invalid argument: 'default_friction_coefficients' must be positive!");
  }
}


inline ContactStatus::ContactStatus() 
  : contact_types_(),
    contact_frame_names_(),
    is_contact_active_(),
    contact_placements_(),
    contact_positions_(),
    contact_rotations_(),
    friction_coefficients_(),
    dimf_(0),
    max_contacts_(0),
    max_num_contacts_(0),
    has_active_contacts_(false) {
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


inline const std::string& ContactStatus::contactFrameName(
    const int contact_index) const {
  assert(contact_index >= 0);
  assert(contact_index < max_num_contacts_);
  if (contact_frame_names_.empty()) {
    throw std::runtime_error("[ContactStatus] invalid argument: contact_frame_names_ is empty!");
  }
  return contact_frame_names_[contact_index];
}


inline const std::vector<std::string>& ContactStatus::contactFrameNames() const {
  return contact_frame_names_;
}


inline bool ContactStatus::isContactActive(const int contact_index) const {
  assert(!is_contact_active_.empty());
  assert(contact_index >= 0);
  assert(contact_index < is_contact_active_.size());
  return is_contact_active_[contact_index];
}


inline bool ContactStatus::isContactActive(
    const std::string& contact_frame_name) const {
  return isContactActive(findContactIndex(contact_frame_name));
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
        break;
      default:
        break;
    }
  }
  setHasActiveContacts();
}


inline void ContactStatus::activateContact(const std::string& contact_frame_name) {
  activateContact(findContactIndex(contact_frame_name));
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
        break;
      default:
        break;
    }
  }
  setHasActiveContacts();
}


inline void ContactStatus::deactivateContact(const std::string& contact_frame_name) {
  deactivateContact(findContactIndex(contact_frame_name));
}


inline void ContactStatus::activateContacts(
    const std::vector<int>& contact_indices) {
  assert(contact_indices.size() <= max_num_contacts_);
  for (const int e : contact_indices) {
    activateContact(e);
  }
}
 

inline void ContactStatus::activateContacts(
    const std::vector<std::string>& contact_frame_names) {
  assert(contact_frame_names.size() <= max_num_contacts_);
  for (const auto& e : contact_frame_names) {
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


inline void ContactStatus::deactivateContacts(
    const std::vector<std::string>& contact_frame_names) {
  assert(contact_frame_names.size() <= max_num_contacts_);
  for (const auto& e : contact_frame_names) {
    deactivateContact(e);
  }
}


inline void ContactStatus::setContactPlacement(
    const int contact_index, const Eigen::Vector3d& contact_position) {
  setContactPlacement(contact_index, contact_position, 
                      Eigen::Matrix3d::Identity());
}


inline void ContactStatus::setContactPlacement(
    const std::string& contact_frame_name,
    const Eigen::Vector3d& contact_position) {
  setContactPlacement(findContactIndex(contact_frame_name), contact_position, 
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


inline void ContactStatus::setContactPlacement(
    const std::string& contact_frame_name, 
    const Eigen::Vector3d& contact_position,
    const Eigen::Matrix3d& contact_rotation) {
  setContactPlacement(findContactIndex(contact_frame_name), contact_position,
                      contact_rotation);
}


inline void ContactStatus::setContactPlacement(const int contact_index, 
                                               const SE3& contact_placement) {
  assert(contact_index >= 0);
  assert(contact_index < max_num_contacts_);
  contact_positions_[contact_index] = contact_placement.translation();
  contact_rotations_[contact_index] = contact_placement.rotation();
  contact_placements_[contact_index] = contact_placement;
}


inline void ContactStatus::setContactPlacement(const std::string& contact_frame_name, 
                                               const SE3& contact_placement) {
  setContactPlacement(contact_frame_name, contact_placement);
}


inline void ContactStatus::setContactPlacements(
    const std::vector<Eigen::Vector3d>& contact_positions) {
  assert(contact_positions.size() == max_num_contacts_);
  for (int i=0; i<max_num_contacts_; ++i) {
    setContactPlacement(i, contact_positions[i], Eigen::Matrix3d::Identity());
  }
}


inline void ContactStatus::setContactPlacements(
    const std::unordered_map<std::string, Eigen::Vector3d>& contact_positions) {
  assert(contact_positions.size() == max_num_contacts_);
  for (int i=0; i<max_num_contacts_; ++i) {
    setContactPlacement(i, contact_positions.at(contactFrameName(i)),
                        Eigen::Matrix3d::Identity());
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


inline void ContactStatus::setContactPlacements(
    const std::unordered_map<std::string, Eigen::Vector3d>& contact_positions,
    const std::unordered_map<std::string, Eigen::Matrix3d>& contact_rotations) {
  assert(contact_positions.size() == max_num_contacts_);
  assert(contact_rotations.size() == max_num_contacts_);
  for (int i=0; i<max_num_contacts_; ++i) {
    setContactPlacement(i, contact_positions.at(contactFrameName(i)),
                        contact_rotations.at(contactFrameName(i)));
  }
}


inline void ContactStatus::setContactPlacements(
    const aligned_vector<SE3>& contact_placements) {
  assert(contact_placements.size() == max_num_contacts_);
  for (int i=0; i<max_num_contacts_; ++i) {
    setContactPlacement(i, contact_placements[i]);
  }
}


inline void ContactStatus::setContactPlacements(
    const aligned_unordered_map<std::string, SE3>& contact_placements) {
  assert(contact_placements.size() == max_num_contacts_);
  for (int i=0; i<max_num_contacts_; ++i) {
    setContactPlacement(i, contact_placements.at(contactFrameName(i)));
  }
}


inline const SE3& ContactStatus::contactPlacement(
    const int contact_index) const {
  assert(contact_index >= 0);
  assert(contact_index < max_num_contacts_);
  return contact_placements_[contact_index];
}


inline const SE3& ContactStatus::contactPlacement(
    const std::string& contact_frame_name) const {
  return contactPlacement(findContactIndex(contact_frame_name));
}


inline const Eigen::Vector3d& ContactStatus::contactPosition(
    const int contact_index) const {
  assert(contact_index >= 0);
  assert(contact_index < max_num_contacts_);
  return contact_positions_[contact_index];
}


inline const Eigen::Vector3d& ContactStatus::contactPosition(
    const std::string& contact_frame_name) const {
  return contactPosition(findContactIndex(contact_frame_name));
}


inline const Eigen::Matrix3d& ContactStatus::contactRotation(
    const int contact_index) const {
  assert(contact_index >= 0);
  assert(contact_index < max_num_contacts_);
  return contact_rotations_[contact_index];
}


inline const Eigen::Matrix3d& ContactStatus::contactRotation(
    const std::string& contact_frame_name) const {
  return contactRotation(findContactIndex(contact_frame_name));
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


inline void ContactStatus::setFrictionCoefficient(const int contact_index,  
                                                  const double friction_coefficient) {
  if (friction_coefficient <= 0.0) {
    throw std::invalid_argument("[ContactStatus] invalid argument: 'friction_coefficient' must be positive!");
  }
  assert(contact_index >= 0);
  assert(contact_index < max_num_contacts_);
  friction_coefficients_[contact_index] = friction_coefficient;
}


inline void ContactStatus::setFrictionCoefficient(const std::string& contact_frame_name, 
                                                  const double friction_coefficient) {
  setFrictionCoefficient(findContactIndex(contact_frame_name), friction_coefficient);
}


inline void ContactStatus::setFrictionCoefficients(
      const std::vector<double>& friction_coefficients) {
  assert(friction_coefficients.size() == max_num_contacts_);
  for (int i=0; i<max_num_contacts_; ++i) {
    setFrictionCoefficient(i, friction_coefficients[i]);
  }
}


inline void ContactStatus::setFrictionCoefficients(
      const std::unordered_map<std::string, double>& friction_coefficients) {
  assert(friction_coefficients.size() == max_num_contacts_);
  for (int i=0; i<max_num_contacts_; ++i) {
    setFrictionCoefficient(i, friction_coefficients.at(contactFrameName(i)));
  }
}


inline double ContactStatus::frictionCoefficient(const int contact_index) const {
  assert(contact_index >= 0);
  assert(contact_index < max_num_contacts_);
  return friction_coefficients_[contact_index];
}


inline double ContactStatus::frictionCoefficient(
    const std::string& contact_frame_name) const {
  return frictionCoefficient(findContactIndex(contact_frame_name));
}


inline const std::vector<double>& ContactStatus::frictionCoefficients() const {
  return friction_coefficients_;
}


inline int ContactStatus::findContactIndex(
    const std::string& contact_frame_name) const {
  if (contact_frame_names_.empty()) {
    throw std::runtime_error("[ContactStatus] invalid argument: contact_frame_names_ is empty!");
  }
  for (int i=0; i<contact_frame_names_.size(); ++i) {
    if (contact_frame_names_[i] == contact_frame_name) {
      return i;
    }
  }
  throw std::runtime_error("[ContactStatus] cannot find the input contact_frame_name: '" + contact_frame_name + "'");
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


inline void ContactStatus::setHasActiveContacts() {
  has_active_contacts_ = (dimf_ > 0);
}


} // namespace robotoc

#endif // ROBOTOC_CONTACT_STATUS_HXX_ 