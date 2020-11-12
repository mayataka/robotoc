#ifndef IDOCP_CONTACT_SEQUENCE_HXX
#define IDOCP_CONTACT_SEQUENCE_HXX_

#include "idocp/ocp/contact_sequence.hpp"

#include <stdexcept>
#include <cassert>
#include <cmath>

namespace idocp {

inline ContactSequence::ContactSequence(const Robot& robot, const double T, 
                                        const int N)
  : N_(N),
    T_(T),
    dtau_(T/N),
    contact_sequence_(N, robot.createContactStatus()),
    discrete_event_sequence_(N, DiscreteEvent(robot)),
    num_impulse_stages_(N+1, 0), 
    num_lift_stages_(N+1, 0),
    impulse_stage_(N, -1),
    lift_stage_(N, -1) {
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
  : N_(0),
    T_(0),
    dtau_(0),
    contact_sequence_(),
    discrete_event_sequence_(),
    num_impulse_stages_(), 
    num_lift_stages_(),
    impulse_stage_(), 
    lift_stage_() {
}


inline ContactSequence::~ContactSequence() {
}


inline void ContactSequence::setContactStatusUniformly(
    const ContactStatus& contact_status) {
  for (auto& e : contact_sequence_) {
    e.set(contact_status);
  }
  for (auto& e : discrete_event_sequence_) {
    e.disableDiscreteEvent();
  }
  countAll();
}


inline void ContactSequence::setDiscreteEvent(
    const DiscreteEvent& discrete_event) {
  assert(discrete_event.eventTime() > 0);
  assert(discrete_event.eventTime() < T_);
  assert(discrete_event.existDiscreteEvent());
  const int time_stage = std::floor(discrete_event.eventTime()/dtau_);
  if (discrete_event_sequence_[time_stage].existDiscreteEvent()) {
    discrete_event.act(contact_sequence_[time_stage+1]);
    discrete_event_sequence_[time_stage].setDiscreteEvent(
        contact_sequence_[time_stage], contact_sequence_[time_stage+1]);
  }
  else {
    discrete_event_sequence_[time_stage] = discrete_event;
  }
  for (int i=time_stage+1; i<N_; ++i) {
    discrete_event_sequence_[time_stage].act(contact_sequence_[i]);
    discrete_event_sequence_[i].disableDiscreteEvent();
  }
  countAll();
}


inline void ContactSequence::shiftImpulse(const int impulse_index, 
                                          const double impulse_time) {
  assert(impulse_index >= 0);
  assert(impulse_index < totalNumImpulseStages());
  shiftDiscreteEvent(timeStageBeforeImpulse(impulse_index), impulse_time);
}


inline void ContactSequence::shiftLift(const int lift_index, 
                                       const double lift_time) {
  assert(lift_index >= 0);
  assert(lift_index < totalNumLiftStages());
  shiftDiscreteEvent(timeStageBeforeLift(lift_index), lift_time);
}


inline const ContactStatus& ContactSequence::contactStatus(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return contact_sequence_[time_stage];
}


inline const ImpulseStatus& ContactSequence::impulseStatus(
    const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < totalNumImpulseStages());
  return discrete_event_sequence_[timeStageBeforeImpulse(impulse_index)].impulseStatus();
}


inline double ContactSequence::impulseTime(const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < totalNumImpulseStages());
  return discrete_event_sequence_[timeStageBeforeImpulse(impulse_index)].eventTime();
}


inline double ContactSequence::liftTime(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < totalNumLiftStages());
  return discrete_event_sequence_[timeStageBeforeLift(lift_index)].eventTime();
}


inline int ContactSequence::numImpulseStages(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return num_impulse_stages_[time_stage];
}


inline int ContactSequence::numLiftStages(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return num_lift_stages_[time_stage];
}


inline int ContactSequence::totalNumImpulseStages() const {
  return num_impulse_stages_[N_];
}


inline int ContactSequence::totalNumLiftStages() const {
  return num_lift_stages_[N_];
}


inline int ContactSequence::timeStageBeforeImpulse(
    const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < totalNumImpulseStages());
  return impulse_stage_[impulse_index];
}


inline int ContactSequence::timeStageBeforeLift(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < totalNumLiftStages());
  return lift_stage_[lift_index];
}


inline void ContactSequence::shiftDiscreteEvent(const int time_stage, 
                                                const double event_time) {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  assert(existDiscreteEvent(time_stage));
  if (event_time < 0) {
    shiftDiscreteEventBeyondInitial(time_stage);
  }
  else if (event_time > T_) {
    shiftDiscreteEventBeyondTerminal(time_stage);
  }
  else {
    const int shifted_time_stage = timeStageFromContinuousTime(event_time);
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
    countAll();
  }
}


inline void ContactSequence::shiftDiscreteEventBeyondInitial(
    const int time_stage) {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  assert(existDiscreteEvent(time_stage));
  for (int i=time_stage; i>=0; --i) {
    contact_sequence_[i].set(contact_sequence_[time_stage+1]);
    discrete_event_sequence_[i].disableDiscreteEvent();
  }
  countAll();
}


inline void ContactSequence::shiftDiscreteEventBeyondTerminal(
    const int time_stage) {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  assert(existDiscreteEvent(time_stage));
  discrete_event_sequence_[time_stage].disableDiscreteEvent();
  for (int i=time_stage+1; i<N_; ++i) {
    contact_sequence_[i].set(contact_sequence_[time_stage]);
    discrete_event_sequence_[i].disableDiscreteEvent();
  }
  countAll();
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
  return (discrete_event_sequence_[time_stage].existLift()
            && (!existImpulse(time_stage)));
}


inline int ContactSequence::timeStageFromContinuousTime(
    const double time) const {
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


inline void ContactSequence::countAll() {
  countNumImpulseStaeges();
  countNumLiftStages();
  countImpulseStage();
  countLiftStage();
}


inline void ContactSequence::countNumImpulseStaeges() {
  num_impulse_stages_[0] = 0;
  for (int i=0; i<N_; ++i) {
    if (existImpulse(i)) num_impulse_stages_[i+1] = num_impulse_stages_[i] + 1;
    else num_impulse_stages_[i+1] = num_impulse_stages_[i];
  }
}


inline void ContactSequence::countNumLiftStages() {
  num_lift_stages_[0] = 0;
  for (int i=0; i<N_; ++i) {
    if (existLift(i)) num_lift_stages_[i+1] = num_lift_stages_[i] + 1;
    else num_lift_stages_[i+1] = num_lift_stages_[i];
  }
}


inline void ContactSequence::countImpulseStage() {
  int iterator = 0;
  for (int i=0; i<N_; ++i) {
    if (existImpulse(i)) {
      impulse_stage_[iterator] = i;
      ++iterator;
    }
  }
  for (int i=iterator; i<N_; ++i) {
    impulse_stage_[i] = -1;
  }
}


inline void ContactSequence::countLiftStage() {
  int iterator = 0;
  for (int i=0; i<N_; ++i) {
    if (existLift(i)) {
      lift_stage_[iterator] = i;
      ++iterator;
    }
  }
  for (int i=iterator; i<N_; ++i) {
    lift_stage_[i] = -1;
  }
}

} // namespace idocp 

#endif // IDOCP_CONTACT_SEQUENCE_HXX_ 