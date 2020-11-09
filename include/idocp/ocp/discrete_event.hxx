#ifndef IDOCP_DISCRETE_EVENT_HXX_
#define IDOCP_DISCRETE_EVENT_HXX_

#include "idocp/ocp/discrete_event.hpp"

#include <cassert>


namespace idocp {

inline DiscreteEvent::DiscreteEvent(const int max_point_contacts)
  : contact_status_before_(max_point_contacts),
    contact_status_after_(max_point_contacts),
    impulse_status_(max_point_contacts),
    max_point_contacts_(max_point_contacts),
    has_impulse_(false), 
    has_lift_(false),
    time_(0) {
}


inline DiscreteEvent::DiscreteEvent() 
  : contact_status_before_(),
    contact_status_after_(),
    impulse_status_(),
    max_point_contacts_(0),
    has_impulse_(false), 
    has_lift_(false),
    time_(0) {
}
 

inline DiscreteEvent::~DiscreteEvent() {
}


inline const ImpulseStatus& DiscreteEvent::impulseStatus() const {
  return impulse_status_;
}


inline bool DiscreteEvent::hasDiscreteEvent() const {
  return (has_impulse_ || has_lift_);
}


inline bool DiscreteEvent::hasImpulse() const {
  return has_impulse_;
}


inline bool DiscreteEvent::hasLift() const {
  return has_lift_;
}


inline double DiscreteEvent::eventTime() const {
  return time_;
}


inline void DiscreteEvent::act(ContactStatus& contact_status) const {
  assert(contact_status.max_point_contacts() == max_point_contacts_);
  assert(contact_status == contact_status_before_);
  contact_status.setContactStatus(contact_status_after_.isContactActive());
}


inline void DiscreteEvent::setDiscreteEvent(
    const ContactStatus& contact_status_before, 
    const ContactStatus& contact_status_after) {
  assert(contact_status_before.max_point_contacts() == max_point_contacts_);
  assert(contact_status_after.max_point_contacts() == max_point_contacts_);
  assert(contact_status_before != contact_status_after);
  has_impulse_ = false;
  has_lift_ = false;
  for (int i=0; i<max_point_contacts_; ++i) {
    if (contact_status_before.isContactActive(i)) {
      impulse_status_.deactivateImpulse(i);
      if (!contact_status_after.isContactActive(i)) {
        has_lift_ = true;
      }
    }
    else {
      if (contact_status_after.isContactActive(i)) {
        impulse_status_.activateImpulse(i);
        has_impulse_ = true;
      }
      else {
        impulse_status_.deactivateImpulse(i);
      }
    }
  }
  contact_status_before_.setContactStatus(contact_status_before.isContactActive());
  contact_status_after_.setContactStatus(contact_status_after.isContactActive());
}


inline void DiscreteEvent::setEventTime(const double time) {
  time_ = time;
}


inline void DiscreteEvent::disableDiscreteEvent() {
  for (int i=0; i<max_point_contacts_; ++i) {
    impulse_status_.deactivateImpulse(i);
  }
  has_impulse_ = false;
  has_lift_ = false;
  time_ = 0;
}


inline int DiscreteEvent::max_point_contacts() const {
  return max_point_contacts_;
}

} // namespace idocp

#endif // IDOCP_DISCRETE_EVENT_HXX_ 