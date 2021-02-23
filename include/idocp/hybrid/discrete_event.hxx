#ifndef IDOCP_DISCRETE_EVENT_HXX_
#define IDOCP_DISCRETE_EVENT_HXX_

#include "idocp/hybrid/discrete_event.hpp"

#include <cassert>


namespace idocp {

inline DiscreteEvent::DiscreteEvent(const int max_point_contacts)
  : pre_contact_status_(max_point_contacts),
    post_contact_status_(max_point_contacts),
    impulse_status_(max_point_contacts),
    max_point_contacts_(max_point_contacts),
    exist_impulse_(false), 
    exist_lift_(false) {
}


inline DiscreteEvent::DiscreteEvent(const ContactStatus& pre_contact_status, 
                                    const ContactStatus& post_contact_status)
  : pre_contact_status_(pre_contact_status.maxPointContacts()),
    post_contact_status_(pre_contact_status.maxPointContacts()),
    impulse_status_(pre_contact_status.maxPointContacts()),
    max_point_contacts_(pre_contact_status.maxPointContacts()),
    exist_impulse_(false), 
    exist_lift_(false) {
  setDiscreteEvent(pre_contact_status, post_contact_status);
}


inline DiscreteEvent::DiscreteEvent() 
  : pre_contact_status_(),
    post_contact_status_(),
    impulse_status_(),
    max_point_contacts_(0),
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
  assert(pre_contact_status.maxPointContacts() == max_point_contacts_);
  assert(post_contact_status.maxPointContacts() == max_point_contacts_);
  exist_impulse_ = false;
  exist_lift_ = false;
  for (int i=0; i<max_point_contacts_; ++i) {
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
  setContactPoints(post_contact_status.contactPoints());
  pre_contact_status_ = pre_contact_status;
  post_contact_status_ = post_contact_status;
}


inline void DiscreteEvent::setContactPoint(
    const int contact_index, const Eigen::Vector3d& contact_point) {
  assert(contact_index >= 0);
  assert(contact_index < max_point_contacts_);
  impulse_status_.setContactPoint(contact_index, contact_point);
}


inline void DiscreteEvent::setContactPoints(
    const std::vector<Eigen::Vector3d>& contact_points) {
  assert(contact_points.size() == max_point_contacts_);
  impulse_status_.setContactPoints(contact_points);
}


inline int DiscreteEvent::maxPointContacts() const {
  return max_point_contacts_;
}

} // namespace idocp

#endif // IDOCP_DISCRETE_EVENT_HXX_ 