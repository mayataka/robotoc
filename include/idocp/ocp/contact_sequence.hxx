#ifndef IDOCP_CONTACT_SEQUENCE_HXX
#define IDOCP_CONTACT_SEQUENCE_HXX_

#include "idocp/ocp/contact_sequence.hpp"

#include <stdexcept>
#include <cassert>
#include <cmath>

namespace idocp {

inline ContactSequence::ContactSequence(const Robot& robot, const double T, 
                                        const int N)
  : max_point_contacts_(robot.max_point_contacts()),
    N_(N),
    T_(T),
    dtau_(T/N),
    contact_sequence_(N, ContactStatus(robot.max_point_contacts())),
    discrete_event_sequence_(N, DiscreteEvent(robot.max_point_contacts())) {
  try {
    if (T <= 0) {
      throw std::out_of_range("invalid value: T must be positive!");
    }
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline ContactSequence::ContactSequence()
  : max_point_contacts_(0),
    N_(0),
    T_(0),
    dtau_(0),
    contact_sequence_(),
    discrete_event_sequence_() {
}


inline ContactSequence::~ContactSequence() {
}


inline void ContactSequence::setContactStatusUniformly(
    const ContactStatus& contact_status) {
  assert(contact_status.max_point_contacts() == max_point_contacts_);
  for (auto e : contact_sequence_) {
    e.isContactActive() = contact_status.isContactActive();
  }
}


inline void ContactSequence::setDiscreteEvent(
    const DiscreteEvent& discrete_event) {
  assert(discrete_event.eventTime() > 0);
  assert(discrete_event.eventTime() < T_);
  assert(discrete_event.hasDiscreteEvent());
  const int time_stage = std::floor(discrete_event.eventTime()/dtau_);
  setDiscreteEvent(time_stage, discrete_event);
}


inline void ContactSequence::setDiscreteEvent(
    const int time_stage, const DiscreteEvent& discrete_event) {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  assert(discrete_event.hasDiscreteEvent());
  discrete_event_sequence_[time_stage] = discrete_event;
  for (int i=time_stage+1; i<N_; ++i) {
    discrete_event.act(contact_sequence_[i]);
    discrete_event_sequence_[i+1].disableDiscreteEvent();
  }
}


inline void ContactSequence::shiftDiscreteEvent(const int time_stage, 
                                                const double shift_event_time) {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  discrete_event_sequence_[time_stage] = discrete_event;
}


inline void ContactSequence::shiftDiscreteEvent(const int time_stage, 
                                                const int shift_time_stages) {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  discrete_event_sequence_[time_stage] = discrete_event;
}


inline const ContactStatus& ContactSequence::contactStatus(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return contact_sequence_[time_stage];
}


inline const ImpulseStatus& ContactSequence::impulseStatus(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return discrete_event_sequence_[time_stage].impulseStatus();
}


inline bool ContactSequence::hasDiscreteEvent(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return (hasImpulse(time_stage) || hasLift(time_stage));
}


inline bool ContactSequence::hasImpulse(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return discrete_event_sequence_[time_stage].hasImpulse();
}


inline bool ContactSequence::hasLift(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return discrete_event_sequence_[time_stage].hasLift();
}


inline void ContactSequence::setContactSequenceFromDiscreteEvent(
    const int time_stage_begin) {
  for (int i=time_stage_begin; i<N_; ++i) {

  }
}

} // namespace idocp 

#endif // IDOCP_CONTACT_SEQUENCE_HXX_ 