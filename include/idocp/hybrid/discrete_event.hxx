#ifndef IDOCP_DISCRETE_EVENT_HXX_
#define IDOCP_DISCRETE_EVENT_HXX_

#include "idocp/hybrid/discrete_event.hpp"

#include <cassert>


namespace idocp {

inline DiscreteEvent::DiscreteEvent(const int max_point_contacts)
  : eventTime(0),
    pre_contact_status_(max_point_contacts),
    post_contact_status_(max_point_contacts),
    impulse_status_(max_point_contacts),
    max_point_contacts_(max_point_contacts),
    exist_impulse_(false), 
    exist_lift_(false) {
}


inline DiscreteEvent::DiscreteEvent(const Robot& robot)
  : eventTime(0),
    pre_contact_status_(robot.max_point_contacts()),
    post_contact_status_(robot.max_point_contacts()),
    impulse_status_(robot.max_point_contacts()),
    max_point_contacts_(robot.max_point_contacts()),
    exist_impulse_(false), 
    exist_lift_(false) {
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


inline const ImpulseStatus& DiscreteEvent::impulseStatus() const {
  return impulse_status_;
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


inline const ContactStatus& DiscreteEvent::preContactStatus() const {
  return pre_contact_status_;
}


inline const ContactStatus& DiscreteEvent::postContactStatus() const {
  return post_contact_status_;
}


inline void DiscreteEvent::setDiscreteEvent(
    const ContactStatus& pre_contact_status, 
    const ContactStatus& post_contact_status) {
  assert(pre_contact_status.max_point_contacts() == max_point_contacts_);
  assert(post_contact_status.max_point_contacts() == max_point_contacts_);
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
  pre_contact_status_.set(pre_contact_status);
  post_contact_status_.set(post_contact_status);
}


inline void DiscreteEvent::setContactPoints(const Robot& robot) {
  robot.setContactPoints(impulse_status_);
}


inline void DiscreteEvent::disableDiscreteEvent() {
  for (int i=0; i<max_point_contacts_; ++i) {
    impulse_status_.deactivateImpulse(i);
    impulse_status_.setContactPoint(i, Eigen::Vector3d::Zero());
  }
  exist_impulse_ = false;
  exist_lift_ = false;
  eventTime = 0;
  post_contact_status_.set(pre_contact_status_);
}


inline int DiscreteEvent::max_point_contacts() const {
  return max_point_contacts_;
}

} // namespace idocp

#endif // IDOCP_DISCRETE_EVENT_HXX_ 