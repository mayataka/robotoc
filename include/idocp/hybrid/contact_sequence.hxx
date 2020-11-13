#ifndef IDOCP_CONTACT_SEQUENCE_HXX
#define IDOCP_CONTACT_SEQUENCE_HXX_

#include "idocp/hybrid/contact_sequence.hpp"

#include <stdexcept>
#include <cassert>
#include <cmath>

namespace idocp {

inline ContactSequence::ContactSequence(const Robot& robot, const double T, 
                                        const int N)
  : N_(N),
    T_(T),
    dtau_(T/N),
    contact_sequence_(robot, N),
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
    contact_sequence_(),
    num_impulse_stages_(), 
    num_lift_stages_(),
    impulse_stage_(), 
    lift_stage_() {
}


inline ContactSequence::~ContactSequence() {
}


inline void ContactSequence::setContactStatusUniformly(
    const ContactStatus& contact_status) {
  contact_sequence_.setContactStatusUniformly(contact_status);
  countAll();
}


inline void ContactSequence::setDiscreteEvent(
    const DiscreteEvent& discrete_event) {
  assert(discrete_event.eventTime > 0);
  assert(discrete_event.eventTime < T_);
  assert(discrete_event.existDiscreteEvent());
  const int event_time_stage = std::floor(discrete_event.eventTime/dtau_);
  contact_sequence_.setDiscreteEvent(discrete_event, event_time_stage);
  contact_sequence_.setEventTime(event_time_stage, discrete_event.eventTime);
  countAll();
}


inline void ContactSequence::shiftImpulse(const int impulse_index, 
                                          const double impulse_time) {
  assert(impulse_index >= 0);
  assert(impulse_index < totalNumImpulseStages());
  const int shifted_impulse_time_stage = std::floor(impulse_time/dtau_);
  contact_sequence_.shiftDiscreteEvent(timeStageBeforeImpulse(impulse_index), 
                                       shifted_impulse_time_stage);
  countAll();
}


inline void ContactSequence::shiftLift(const int lift_index, 
                                       const double lift_time) {
  assert(lift_index >= 0);
  assert(lift_index < totalNumImpulseStages());
  const int shifted_lift_time_stage = std::floor(lift_time/dtau_);
  contact_sequence_.shiftDiscreteEvent(timeStageBeforeLift(lift_index), 
                                       shifted_lift_time_stage);
  countAll();
}


inline const ContactStatus& ContactSequence::contactStatus(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return contact_sequence_.contactStatus(time_stage);
}


inline const ImpulseStatus& ContactSequence::impulseStatus(
    const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < totalNumImpulseStages());
  return contact_sequence_.impulseStatus(timeStageBeforeImpulse(impulse_index));
}


inline double ContactSequence::impulseTime(const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < totalNumImpulseStages());
  return contact_sequence_.eventTime(timeStageBeforeImpulse(impulse_index));
}


inline double ContactSequence::liftTime(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < totalNumLiftStages());
  return contact_sequence_.eventTime(timeStageBeforeLift(lift_index));
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
    if (contact_sequence_.existImpulse(i)) 
      num_impulse_stages_[i+1] = num_impulse_stages_[i] + 1;
    else 
      num_impulse_stages_[i+1] = num_impulse_stages_[i];
  }
}


inline void ContactSequence::countNumLiftStages() {
  num_lift_stages_[0] = 0;
  for (int i=0; i<N_; ++i) {
    if (contact_sequence_.existOnlyLift(i)) 
      num_lift_stages_[i+1] = num_lift_stages_[i] + 1;
    else 
      num_lift_stages_[i+1] = num_lift_stages_[i];
  }
}


inline void ContactSequence::countImpulseStage() {
  int iterator = 0;
  for (int i=0; i<N_; ++i) {
    if (contact_sequence_.existImpulse(i)) {
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
    if (contact_sequence_.existOnlyLift(i)) {
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