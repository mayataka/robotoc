#ifndef ROBOTOC_DISCRETE_EVENT_HXX_
#define ROBOTOC_DISCRETE_EVENT_HXX_

#include "robotoc/hybrid/discrete_event.hpp"

#include <cassert>


namespace robotoc {

inline DiscreteEvent::DiscreteEvent(
    const std::vector<ContactType>& contact_types)
  : pre_contact_status_(contact_types),
    post_contact_status_(contact_types),
    impulse_status_(contact_types),
    max_num_contacts_(contact_types.size()),
    event_type_(DiscreteEventType::None),
    exist_impulse_(false), 
    exist_lift_(false) {
}


inline DiscreteEvent::DiscreteEvent(const ContactStatus& pre_contact_status, 
                                    const ContactStatus& post_contact_status)
  : pre_contact_status_(pre_contact_status.contactTypes()),
    post_contact_status_(pre_contact_status.contactTypes()),
    impulse_status_(pre_contact_status.contactTypes()),
    max_num_contacts_(pre_contact_status.maxNumContacts()),
    event_type_(DiscreteEventType::None),
    exist_impulse_(false), 
    exist_lift_(false) {
  setDiscreteEvent(pre_contact_status, post_contact_status);
}


inline DiscreteEvent::DiscreteEvent() 
  : pre_contact_status_(),
    post_contact_status_(),
    impulse_status_(),
    max_num_contacts_(0),
    event_type_(DiscreteEventType::None),
    exist_impulse_(false), 
    exist_lift_(false) {
}
 

inline DiscreteEvent::~DiscreteEvent() {
}


inline bool DiscreteEvent::existDiscreteEvent() const {
  return (exist_impulse_ || exist_lift_);
}


inline bool DiscreteEvent::existImpulse() const {
  return exist_impulse_;
}


inline bool DiscreteEvent::existLift() const {
  return exist_lift_;
}


inline const ImpulseStatus& DiscreteEvent::impulseStatus() const {
  return impulse_status_;
}


inline const ContactStatus& DiscreteEvent::preContactStatus() const {
  return pre_contact_status_;
}


inline const ContactStatus& DiscreteEvent::postContactStatus() const {
  return post_contact_status_;
}


inline void DiscreteEvent::setDiscreteEvent(
    const ContactStatus& pre_contact_status, 
    const ContactStatus& post_contact_status) {
  assert(pre_contact_status.maxNumContacts() == max_num_contacts_);
  assert(post_contact_status.maxNumContacts() == max_num_contacts_);
  exist_impulse_ = false;
  exist_lift_ = false;
  for (int i=0; i<max_num_contacts_; ++i) {
    if (pre_contact_status.isContactActive(i)) {
      impulse_status_.deactivateImpulse(i);
      if (!post_contact_status.isContactActive(i)) {
        exist_lift_ = true;
      }
    }
    else {
      if (post_contact_status.isContactActive(i)) {
        impulse_status_.activateImpulse(i);
        exist_impulse_ = true;
      }
      else {
        impulse_status_.deactivateImpulse(i);
      }
    }
  }
  impulse_status_.setImpulseModeId(pre_contact_status.contactModeId());
  setContactPlacements(post_contact_status.contactPositions(),
                       post_contact_status.contactRotations());
  pre_contact_status_ = pre_contact_status;
  post_contact_status_ = post_contact_status;
  if (exist_impulse_) { event_type_ = DiscreteEventType::Impulse; }
  else if (exist_lift_) { event_type_ = DiscreteEventType::Lift; }
  else { event_type_ = DiscreteEventType::None; }
}


inline void DiscreteEvent::setContactPlacement(
    const int contact_index, const Eigen::Vector3d& contact_position) {
  assert(contact_index >= 0);
  assert(contact_index < max_num_contacts_);
  impulse_status_.setContactPlacement(contact_index, contact_position);
}


inline void DiscreteEvent::setContactPlacement(
    const int contact_index, const Eigen::Vector3d& contact_position,
    const Eigen::Matrix3d& contact_rotation) {
  assert(contact_index >= 0);
  assert(contact_index < max_num_contacts_);
  impulse_status_.setContactPlacement(contact_index, contact_position, 
                                      contact_rotation);
}


inline void DiscreteEvent::setContactPlacements(
    const std::vector<Eigen::Vector3d>& contact_positions) {
  impulse_status_.setContactPlacements(contact_positions);
}


inline void DiscreteEvent::setContactPlacements(
    const std::vector<Eigen::Vector3d>& contact_positions,
    const std::vector<Eigen::Matrix3d>& contact_rotations) {
  impulse_status_.setContactPlacements(contact_positions, contact_rotations);
}


inline void DiscreteEvent::setContactPlacements(
    const aligned_vector<SE3>& contact_placements) {
  impulse_status_.setContactPlacements(contact_placements);
}


inline int DiscreteEvent::maxNumContacts() const {
  return max_num_contacts_;
}


inline DiscreteEventType DiscreteEvent::eventType() const {
  return event_type_;
}

} // namespace robotoc

#endif // ROBOTOC_DISCRETE_EVENT_HXX_ 