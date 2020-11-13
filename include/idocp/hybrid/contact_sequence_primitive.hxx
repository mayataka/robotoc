#ifndef IDOCP_CONTACT_SEQUENCE_PRIMITIVE_HXX_
#define IDOCP_CONTACT_SEQUENCE_PRIMITIVE_HXX_ 

#include "idocp/hybrid/contact_sequence_primitive.hpp"

#include <stdexcept>
#include <cassert>
#include <cmath>

namespace idocp {

inline ContactSequencePrimitive::ContactSequencePrimitive(const Robot& robot, 
                                                          const int N)
  : N_(N),
    contact_sequence_(N+1, robot.createContactStatus()),
    discrete_event_sequence_(N, DiscreteEvent(robot)) {
  try {
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline ContactSequencePrimitive::ContactSequencePrimitive()
  : N_(0),
    contact_sequence_(),
    discrete_event_sequence_() {
}


inline ContactSequencePrimitive::~ContactSequencePrimitive() {
}


inline void ContactSequencePrimitive::setContactStatusUniformly(
    const ContactStatus& contact_status) {
  for (auto& e : contact_sequence_) {
    e.set(contact_status);
  }
  for (auto& e : discrete_event_sequence_) {
    e.disableDiscreteEvent();
  }
}


inline void ContactSequencePrimitive::setDiscreteEvent(
    const DiscreteEvent& discrete_event, const int time_stage) {
  assert(discrete_event.existDiscreteEvent());
  assert(discrete_event.preContactStatus() == contact_sequence_[time_stage]);
  if (discrete_event_sequence_[time_stage].existDiscreteEvent()) {
    contact_sequence_[time_stage+1].set(discrete_event.postContactStatus());
    discrete_event_sequence_[time_stage].setDiscreteEvent(
        contact_sequence_[time_stage], contact_sequence_[time_stage+1]);
  }
  else {
    discrete_event_sequence_[time_stage] = discrete_event;
  }
  for (int i=time_stage+1; i<=N_; ++i) {
    contact_sequence_[i].set(
        discrete_event_sequence_[time_stage].postContactStatus());
  }
  for (int i=time_stage+1; i<N_; ++i) {
    discrete_event_sequence_[i].disableDiscreteEvent();
  }
}


inline void ContactSequencePrimitive::shiftDiscreteEvent(
    const int event_time_stage, const int shifted_event_time_stage) {
  assert(event_time_stage >= 0);
  assert(event_time_stage < N_);
  assert(existDiscreteEvent(event_time_stage));
  if (shifted_event_time_stage < 0) {
    shiftDiscreteEventBeyondInitial(event_time_stage);
  }
  else if (shifted_event_time_stage >= N_) {
    shiftDiscreteEventBeyondTerminal(event_time_stage);
  }
  else {
    if (shifted_event_time_stage < event_time_stage) {
      for (int i=event_time_stage; i>=shifted_event_time_stage+1; --i) {
        contact_sequence_[i].set(contact_sequence_[event_time_stage+1]);
      }
      for (int i=event_time_stage; i>=shifted_event_time_stage+1; --i) {
        discrete_event_sequence_[i].disableDiscreteEvent();
      }
    }
    else if (shifted_event_time_stage > event_time_stage) {
      for (int i=event_time_stage+1; i<=shifted_event_time_stage; ++i) {
        contact_sequence_[i].set(contact_sequence_[event_time_stage]);
      }
      for (int i=event_time_stage; i<shifted_event_time_stage; ++i) {
        discrete_event_sequence_[i].disableDiscreteEvent();
      }
    }
    if (shifted_event_time_stage != event_time_stage) {
      discrete_event_sequence_[shifted_event_time_stage].setDiscreteEvent(
          contact_sequence_[shifted_event_time_stage],
          contact_sequence_[shifted_event_time_stage+1]);
    }
  }
}


inline void ContactSequencePrimitive::shiftDiscreteEventBeyondInitial(
    const int event_time_stage) {
  assert(event_time_stage >= 0);
  assert(event_time_stage < N_);
  assert(existDiscreteEvent(event_time_stage));
  for (int i=event_time_stage; i>=0; --i) {
    contact_sequence_[i].set(contact_sequence_[event_time_stage+1]);
  }
  for (int i=event_time_stage; i>=0; --i) {
    discrete_event_sequence_[i].disableDiscreteEvent();
  }
}


inline void ContactSequencePrimitive::shiftDiscreteEventBeyondTerminal(
    const int event_time_stage) {
  assert(event_time_stage >= 0);
  assert(event_time_stage < N_);
  assert(existDiscreteEvent(event_time_stage));
  for (int i=event_time_stage+1; i<=N_; ++i) {
    contact_sequence_[i].set(contact_sequence_[event_time_stage]);
  }
  for (int i=event_time_stage; i<N_; ++i) {
    discrete_event_sequence_[i].disableDiscreteEvent();
  }
}


inline const ContactStatus& ContactSequencePrimitive::contactStatus(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N_);
  return contact_sequence_[time_stage];
}


inline const ImpulseStatus& ContactSequencePrimitive::impulseStatus(
    const int event_time_stage) const {
  assert(event_time_stage >= 0);
  assert(event_time_stage < N_);
  return discrete_event_sequence_[event_time_stage].impulseStatus();
}


inline double ContactSequencePrimitive::eventTime(
    const int event_time_stage) const {
  assert(event_time_stage >= 0);
  assert(event_time_stage < N_);
  return discrete_event_sequence_[event_time_stage].eventTime;
}


inline void ContactSequencePrimitive::setEventTime(const int event_time_stage, 
                                                   const double event_time) {
  assert(event_time_stage >= 0);
  assert(event_time_stage < N_);
  discrete_event_sequence_[event_time_stage].eventTime = event_time;
}


inline bool ContactSequencePrimitive::existDiscreteEvent(
    const int event_time_stage) const {
  assert(event_time_stage >= 0);
  assert(event_time_stage < N_);
  return (existImpulse(event_time_stage) || existLift(event_time_stage));
}


inline bool ContactSequencePrimitive::existImpulse(
    const int event_time_stage) const {
  assert(event_time_stage >= 0);
  assert(event_time_stage < N_);
  return discrete_event_sequence_[event_time_stage].existImpulse();
}


inline bool ContactSequencePrimitive::existLift(
    const int event_time_stage) const {
  assert(event_time_stage >= 0);
  assert(event_time_stage < N_);
  return discrete_event_sequence_[event_time_stage].existLift();
}


inline bool ContactSequencePrimitive::existOnlyLift(
    const int event_time_stage) const {
  assert(event_time_stage >= 0);
  assert(event_time_stage < N_);
  return ((existLift(event_time_stage)) && (!existImpulse(event_time_stage)));
}

} // namespace idocp 

#endif // IDOCP_CONTACT_SEQUENCE_PRIMITIVE_HXX_ 