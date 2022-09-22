#ifndef ROBOTOC_IMPULSE_STATUS_HXX_
#define ROBOTOC_IMPULSE_STATUS_HXX_

#include "robotoc/robot/impulse_status.hpp"

#include <cassert>
#include <random>


namespace robotoc {

inline ImpulseStatus::ImpulseStatus(
    const std::vector<ContactType>& contact_types, 
    const std::vector<std::string>& contact_frame_names,
    const double default_friction_coefficients,
    const int impulse_mode_id)
  : contact_status_(contact_types, contact_frame_names, 
                    default_friction_coefficients, impulse_mode_id) {
}


inline ImpulseStatus::ImpulseStatus() 
  : contact_status_() {
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


inline int ImpulseStatus::dimf() const {
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


inline void ImpulseStatus::setFrictionCoefficient(const int contact_index, 
                                                  const double friction_coefficient) {
  contact_status_.setFrictionCoefficient(contact_index, friction_coefficient);
}


inline void ImpulseStatus::setFrictionCoefficient(const std::string& contact_frame_name, 
                                                  const double friction_coefficient) {
  contact_status_.setFrictionCoefficient(contact_frame_name, friction_coefficient);
}


inline void ImpulseStatus::setFrictionCoefficients(
    const std::vector<double>& friction_coefficients) {
  contact_status_.setFrictionCoefficients(friction_coefficients);
}


inline void ImpulseStatus::setFrictionCoefficients(
    const std::unordered_map<std::string, double>& friction_coefficients) {
  contact_status_.setFrictionCoefficients(friction_coefficients);
}


inline double ImpulseStatus::frictionCoefficient(const int contact_index) const {
  return contact_status_.frictionCoefficient(contact_index);
}


inline double ImpulseStatus::frictionCoefficient(
    const std::string& contact_frame_name) const {
  return contact_status_.frictionCoefficient(contact_frame_name);
}


inline const std::vector<double>& ImpulseStatus::frictionCoefficients() const {
  return contact_status_.frictionCoefficients();
}


inline int ImpulseStatus::impulseModeId() const {
  return contact_status_.contactModeId();
}


inline void ImpulseStatus::setRandom() {
  contact_status_.setRandom();
}

} // namespace robotoc

#endif // ROBOTOC_IMPULSE_STATUS_HXX_ 