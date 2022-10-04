#ifndef ROBOTOC_IMPACT_STATUS_HXX_
#define ROBOTOC_IMPACT_STATUS_HXX_

#include "robotoc/robot/impact_status.hpp"

#include <cassert>
#include <random>


namespace robotoc {

inline ImpactStatus::ImpactStatus(
    const std::vector<ContactType>& contact_types, 
    const std::vector<std::string>& contact_frame_names,
    const double default_friction_coefficients)
  : contact_status_(contact_types, contact_frame_names, 
                    default_friction_coefficients) {
}


inline ImpactStatus::ImpactStatus() 
  : contact_status_() {
}
 

inline bool ImpactStatus::operator==(const ImpactStatus& other) const {
  assert(other.maxNumContacts() == maxNumContacts());
  for (int i=0; i<maxNumContacts(); ++i) {
    if (other.isImpactActive(i) != isImpactActive(i)) {
      return false;
    }
    if (!other.contactPlacement(i).isApprox(contactPlacement(i))) {
      return false;
    }
  }
  return true;
}


inline bool ImpactStatus::operator!=(const ImpactStatus& other) const {
  return !(*this == other);
}


inline ContactType ImpactStatus::contactType(const int contact_index) const {
  return contact_status_.contactType(contact_index);
}


inline const std::vector<ContactType>& ImpactStatus::contactTypes() const {
  return contact_status_.contactTypes();
}


inline bool ImpactStatus::isImpactActive(const int contact_index) const {
  return contact_status_.isContactActive(contact_index);
}


inline const std::vector<bool>& ImpactStatus::isImpactActive() const {
  return contact_status_.isContactActive();
}


inline bool ImpactStatus::hasActiveImpact() const {
  return contact_status_.hasActiveContacts();
}


inline int ImpactStatus::dimf() const {
  return contact_status_.dimf();
}


inline int ImpactStatus::maxNumContacts() const {
  return contact_status_.maxNumContacts();
}


inline void ImpactStatus::activateImpact(const int impact_index) {
  contact_status_.activateContact(impact_index);
}


inline void ImpactStatus::deactivateImpact(const int impact_index) {
  contact_status_.deactivateContact(impact_index);
}


inline void ImpactStatus::activateImpacts(
    const std::vector<int>& impact_indices) {
  contact_status_.activateContacts(impact_indices);
} 


inline void ImpactStatus::deactivateImpacts(
    const std::vector<int>& impact_indices) {
  contact_status_.deactivateContacts(impact_indices);
}


inline void ImpactStatus::setContactPlacement(
    const int contact_index, const Eigen::Vector3d& contact_position) {
  contact_status_.setContactPlacement(contact_index, contact_position);
}


inline void ImpactStatus::setContactPlacement(
    const int contact_index, const Eigen::Vector3d& contact_position,
    const Eigen::Matrix3d& contact_rotation) {
  contact_status_.setContactPlacement(contact_index, contact_position, 
                                      contact_rotation);
}


inline void ImpactStatus::setContactPlacement(const int contact_index, 
                                               const SE3& contact_placement) {
  contact_status_.setContactPlacement(contact_index, contact_placement);
}


inline void ImpactStatus::setContactPlacements(
    const std::vector<Eigen::Vector3d>& contact_positions) {
  contact_status_.setContactPlacements(contact_positions);
}


inline void ImpactStatus::setContactPlacements(
    const std::vector<Eigen::Vector3d>& contact_positions,
    const std::vector<Eigen::Matrix3d>& contact_rotations) {
  contact_status_.setContactPlacements(contact_positions, contact_rotations);
}


inline void ImpactStatus::setContactPlacements(
    const aligned_vector<SE3>& contact_placements) {
  contact_status_.setContactPlacements(contact_placements);
}


inline const SE3& ImpactStatus::contactPlacement(
    const int contact_index) const {
  return contact_status_.contactPlacement(contact_index);
}


inline const Eigen::Vector3d& ImpactStatus::contactPosition(
    const int contact_index) const {
  return contact_status_.contactPosition(contact_index);
}


inline const Eigen::Matrix3d& ImpactStatus::contactRotation(
    const int contact_index) const {
  return contact_status_.contactRotation(contact_index);
}


inline const aligned_vector<SE3>& ImpactStatus::contactPlacements() const {
  return contact_status_.contactPlacements();
}


inline const std::vector<Eigen::Vector3d>& 
ImpactStatus::contactPositions() const {
  return contact_status_.contactPositions();
}

 
inline const std::vector<Eigen::Matrix3d>& 
ImpactStatus::contactRotations() const {
  return contact_status_.contactRotations();
}


inline void ImpactStatus::setFrictionCoefficient(const int contact_index, 
                                                  const double friction_coefficient) {
  contact_status_.setFrictionCoefficient(contact_index, friction_coefficient);
}


inline void ImpactStatus::setFrictionCoefficient(const std::string& contact_frame_name, 
                                                  const double friction_coefficient) {
  contact_status_.setFrictionCoefficient(contact_frame_name, friction_coefficient);
}


inline void ImpactStatus::setFrictionCoefficients(
    const std::vector<double>& friction_coefficients) {
  contact_status_.setFrictionCoefficients(friction_coefficients);
}


inline void ImpactStatus::setFrictionCoefficients(
    const std::unordered_map<std::string, double>& friction_coefficients) {
  contact_status_.setFrictionCoefficients(friction_coefficients);
}


inline double ImpactStatus::frictionCoefficient(const int contact_index) const {
  return contact_status_.frictionCoefficient(contact_index);
}


inline double ImpactStatus::frictionCoefficient(
    const std::string& contact_frame_name) const {
  return contact_status_.frictionCoefficient(contact_frame_name);
}


inline const std::vector<double>& ImpactStatus::frictionCoefficients() const {
  return contact_status_.frictionCoefficients();
}


inline void ImpactStatus::setRandom() {
  contact_status_.setRandom();
}

} // namespace robotoc

#endif // ROBOTOC_IMPACT_STATUS_HXX_ 