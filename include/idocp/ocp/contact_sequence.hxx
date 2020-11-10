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
  for (auto& e : contact_sequence_) {
    e.set(contact_status);
  }
  for (auto& e : discrete_event_sequence_) {
    e.disableDiscreteEvent();
  }
}


inline void ContactSequence::setDiscreteEvent(
    const DiscreteEvent& discrete_event) {
  assert(discrete_event.eventTime() > 0);
  assert(discrete_event.eventTime() < T_);
  assert(discrete_event.existDiscreteEvent());
  const int time_stage = std::floor(discrete_event.eventTime()/dtau_);
  discrete_event_sequence_[time_stage] = discrete_event;
  for (int i=time_stage+1; i<N_; ++i) {
    discrete_event.act(contact_sequence_[i]);
    discrete_event_sequence_[i].disableDiscreteEvent();
  }
}


inline void ContactSequence::shiftDiscreteEvent(const int time_stage, 
                                                const double event_time) {
  assert(existDiscreteEvent(time_stage));
  if (event_time < 0) {
    shiftDiscreteEventToInitial(time_stage);
  }
  else if (event_time > T_) {
    shiftDiscreteEventToTerminal(time_stage);
  }
  else {
    const int shifted_time_stage = timeStageFromTime(event_time);
    if (shifted_time_stage < time_stage) {
      for (int i=time_stage; i>=shifted_time_stage+1; --i) {
        contact_sequence_[i].set(contact_sequence_[time_stage+1]);
        discrete_event_sequence_[i].disableDiscreteEvent();
      }
    }
    else if (shifted_time_stage > time_stage) {
      for (int i=time_stage; i<=shifted_time_stage; ++i) {
        contact_sequence_[i].set(contact_sequence_[time_stage]);
        discrete_event_sequence_[i].disableDiscreteEvent();
      }
    }
    discrete_event_sequence_[shifted_time_stage].setDiscreteEvent(
        contact_sequence_[shifted_time_stage],
        contact_sequence_[shifted_time_stage+1]);
    if (discrete_event_sequence_[shifted_time_stage].existDiscreteEvent()) {
      discrete_event_sequence_[shifted_time_stage].setEventTime(event_time);
    }
    else {
      discrete_event_sequence_[shifted_time_stage].setEventTime(0);
    }
  }
}


inline void ContactSequence::shiftDiscreteEventToInitial(
    const int time_stage) {
  assert(existDiscreteEvent(time_stage));
  for (int i=time_stage; i>=0; --i) {
    contact_sequence_[i].set(contact_sequence_[time_stage+1]);
    discrete_event_sequence_[i].disableDiscreteEvent();
  }
}


inline void ContactSequence::shiftDiscreteEventToTerminal(
    const int time_stage) {
  assert(existDiscreteEvent(time_stage));
  discrete_event_sequence_[time_stage].disableDiscreteEvent();
  for (int i=time_stage+1; i<N_; ++i) {
    contact_sequence_[i].set(contact_sequence_[time_stage]);
    discrete_event_sequence_[i].disableDiscreteEvent();
  }
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


inline double ContactSequence::eventTime(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return discrete_event_sequence_[time_stage].eventTime();
}


inline bool ContactSequence::existDiscreteEvent(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return (existImpulse(time_stage) || existLift(time_stage));
}


inline bool ContactSequence::existImpulse(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return discrete_event_sequence_[time_stage].existImpulse();
}


inline bool ContactSequence::existLift(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return discrete_event_sequence_[time_stage].existLift();
}


inline int ContactSequence::timeStageFromTime(const double time) const {
  if (time < 0) {
    return 0;
  }
  else if (time > T_) {
    return N_;
  }
  else {
    return std::floor(time/dtau_);
  }
}

} // namespace idocp 

#endif // IDOCP_CONTACT_SEQUENCE_HXX_ 