#ifndef IDOCP_PARNMPC_DISCRETIZER_HXX_
#define IDOCP_PARNMPC_DISCRETIZER_HXX_

#include "idocp/hybrid/parnmpc_discretizer.hpp"

#include <cmath>
#include <cassert>
#include <algorithm>

namespace idocp {

inline ParNMPCDiscretizer::ParNMPCDiscretizer(const double T, const int N, 
                                              const int max_events) 
  : T_(T),
    N_(N),
    N_ideal_(N),
    max_events_(max_events),
    num_impulse_stages_(0),
    num_lift_stages_(0),
    contact_phase_index_from_time_stage_(N+1, 0), 
    impulse_index_before_time_stage_(N, -1), 
    lift_index_before_time_stage_(N, -1), 
    time_stage_after_impulse_(max_events, -1), 
    time_stage_after_lift_(max_events, -1),
    is_time_stage_after_impulse_(N, false),
    is_time_stage_after_lift_(N, false),
    t_(N, 0),
    t_impulse_(max_events, 0),
    t_lift_(max_events, 0),
    dtau_(N, (double)(T/N)),
    dtau_aux_(max_events, 0),
    dtau_lift_(max_events, 0) {
}


inline ParNMPCDiscretizer::ParNMPCDiscretizer()
  : T_(0),
    N_(0),
    N_ideal_(0),
    max_events_(0),
    num_impulse_stages_(0),
    num_lift_stages_(0),
    contact_phase_index_from_time_stage_(), 
    impulse_index_before_time_stage_(), 
    lift_index_before_time_stage_(), 
    time_stage_after_impulse_(), 
    time_stage_after_lift_(),
    is_time_stage_after_impulse_(),
    is_time_stage_after_lift_(),
    t_(),
    t_impulse_(),
    t_lift_(),
    dtau_(),
    dtau_aux_(),
    dtau_lift_() {
}


inline ParNMPCDiscretizer::~ParNMPCDiscretizer() {
}


inline void ParNMPCDiscretizer::discretizeOCP(
    const ContactSequence& contact_sequence, const double t) {
  countImpulseEvents(contact_sequence, t);
  countLiftEvents(contact_sequence, t);
  countIsTimeStageBeforeEvents(contact_sequence);
  countContactPhase(contact_sequence);
  countTime(contact_sequence, t);
}


inline int ParNMPCDiscretizer::N() const {
  return N_;
}


inline int ParNMPCDiscretizer::N_impulse() const {
  return num_impulse_stages_;
}


inline int ParNMPCDiscretizer::N_lift() const {
  return num_lift_stages_;
}


inline int ParNMPCDiscretizer::N_all() const {
  return N()+N_impulse()+N_lift();
}


inline int ParNMPCDiscretizer::N_ideal() const {
  return N_ideal_;
}


inline int ParNMPCDiscretizer::numImpulseStages() const {
  return num_impulse_stages_;
}


inline int ParNMPCDiscretizer::numLiftStages() const {
  return num_lift_stages_;
}


inline int ParNMPCDiscretizer::contactPhase(const int time_stage) const {
  assert(time_stage >= -1);
  assert(time_stage <= N_);
  // When time_stage = -1, i.e., contact_phase_index_from_time_stage_[0] means 
  // the contact state at the initial time of the horizon.
  return contact_phase_index_from_time_stage_[time_stage+1];
}


inline int ParNMPCDiscretizer::impulseIndexBeforeTimeStage(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return impulse_index_before_time_stage_[time_stage];
}


inline int ParNMPCDiscretizer::impulseIndexAfterTimeStage(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return impulseIndexBeforeTimeStage(time_stage+1);
}


inline int ParNMPCDiscretizer::liftIndexBeforeTimeStage(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return lift_index_before_time_stage_[time_stage];
}


inline int ParNMPCDiscretizer::liftIndexAfterTimeStage(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return liftIndexBeforeTimeStage(time_stage+1);
}


inline int ParNMPCDiscretizer::timeStageBeforeImpulse(
    const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < numImpulseStages());
  return timeStageAfterImpulse(impulse_index) - 1;
}


inline int ParNMPCDiscretizer::timeStageAfterImpulse(
    const int impulse_index) const {
  return time_stage_after_impulse_[impulse_index];
}


inline int ParNMPCDiscretizer::timeStageBeforeLift(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < numLiftStages());
  return timeStageAfterLift(lift_index) - 1;
}


inline int ParNMPCDiscretizer::timeStageAfterLift(const int lift_index) const {
  return time_stage_after_lift_[lift_index];
}


inline int ParNMPCDiscretizer::contactPhaseBeforeImpulse(
    const int impulse_index) const {
  return contactPhase(timeStageBeforeImpulse(impulse_index));
}


inline int ParNMPCDiscretizer::contactPhaseAfterImpulse(
    const int impulse_index) const {
  return contactPhase(timeStageAfterImpulse(impulse_index));
}


inline int ParNMPCDiscretizer::contactPhaseBeforeLift(
    const int lift_index) const {
  return contactPhase(timeStageBeforeLift(lift_index));
}


inline int ParNMPCDiscretizer::contactPhaseAfterLift(
    const int lift_index) const {
  return contactPhase(timeStageAfterLift(lift_index));
}


inline bool ParNMPCDiscretizer::isTimeStageBeforeImpulse(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return isTimeStageAfterImpulse(time_stage+1);
}


inline bool ParNMPCDiscretizer::isTimeStageAfterImpulse(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return is_time_stage_after_impulse_[time_stage];
}


inline bool ParNMPCDiscretizer::isTimeStageBeforeLift(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return isTimeStageAfterLift(time_stage+1);
}


inline bool ParNMPCDiscretizer::isTimeStageAfterLift(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return is_time_stage_after_lift_[time_stage];
}


inline double ParNMPCDiscretizer::t(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return t_[time_stage];
}


inline double ParNMPCDiscretizer::t_impulse(const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < numImpulseStages());
  return t_impulse_[impulse_index];
}


inline double ParNMPCDiscretizer::t_lift(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < numLiftStages());
  return t_lift_[lift_index];
}


inline double ParNMPCDiscretizer::dtau(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N_);
  return dtau_[time_stage];
}


inline double ParNMPCDiscretizer::dtau_aux(const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < numImpulseStages());
  return dtau_aux_[impulse_index];
}


inline double ParNMPCDiscretizer::dtau_lift(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < numLiftStages());
  return dtau_lift_[lift_index];
}


inline void ParNMPCDiscretizer::countImpulseEvents(
    const ContactSequence& contact_sequence, const double t) {
  for (auto& e : impulse_index_before_time_stage_) { e = -1; }
  for (auto& e : time_stage_after_impulse_)      { e = -1; }
  num_impulse_stages_ = contact_sequence.numImpulseEvents();
  assert(num_impulse_stages_ <= max_events_);
  const double dt = T_ / N_;
  for (int impulse_index=0; impulse_index<num_impulse_stages_; ++impulse_index) {
    t_impulse_[impulse_index] = contact_sequence.impulseTime(impulse_index) - t;
    time_stage_after_impulse_[impulse_index] 
        = std::floor(t_impulse_[impulse_index]/dt);
    impulse_index_before_time_stage_[time_stage_after_impulse_[impulse_index]] 
        = impulse_index;
  }
}


inline void ParNMPCDiscretizer::countLiftEvents(
    const ContactSequence& contact_sequence, const double t) {
  for (auto& e : lift_index_before_time_stage_) { e = -1; }
  for (auto& e : time_stage_after_lift_)      { e = -1; }
  num_lift_stages_ = contact_sequence.numLiftEvents();
  assert(num_lift_stages_ <=  max_events_);
  const double dt = T_ / N_;
  for (int lift_index=0; lift_index<num_lift_stages_; ++lift_index) {
    t_lift_[lift_index] = contact_sequence.liftTime(lift_index) - t;
    time_stage_after_lift_[lift_index] = std::floor(t_lift_[lift_index]/dt);
    lift_index_before_time_stage_[time_stage_after_lift_[lift_index]] 
        = lift_index;
  }
}


inline void ParNMPCDiscretizer::countIsTimeStageBeforeEvents(
    const ContactSequence& contact_sequence) {
  int impulse_index = 0;
  int lift_index = 0;
  for (int i=0; i<N_; ++i) {
    if (impulse_index < numImpulseStages()) {
      if (i == timeStageAfterImpulse(impulse_index)) {
        is_time_stage_after_impulse_[i] = true;
        ++impulse_index;
      }
      else {
        is_time_stage_after_impulse_[i] = false;
      }
    }
    else {
      is_time_stage_after_impulse_[i] = false;
    }
    if (lift_index < numLiftStages()) {
      if (i == timeStageAfterLift(lift_index)) {
        is_time_stage_after_lift_[i] = true;
        ++lift_index;
      }
      else {
        is_time_stage_after_lift_[i] = false;
      }
    }
    else {
      is_time_stage_after_lift_[i] = false;
    }
  }
}


inline void ParNMPCDiscretizer::countContactPhase(
    const ContactSequence& contact_sequence) {
  for (auto & e : contact_phase_index_from_time_stage_) { e = 0; }
  int num_events = 0;
  for (int i=0; i<N_; ++i) {
    if (isTimeStageAfterImpulse(i) || isTimeStageAfterLift(i)) {
      ++num_events; 
    }
    contact_phase_index_from_time_stage_[i+1] = num_events;
  }
}


inline void ParNMPCDiscretizer::countTime(
    const ContactSequence& contact_sequence, const double t) {
  const double dt = T_ / N_;
  for (int i=0; i<N_; ++i) {
    t_[i] = t + (i+1) * dt;
    dtau_[i] = dt;
  }
  for (int i=0; i<num_impulse_stages_; ++i) {
    dtau_[timeStageAfterImpulse(i)] 
        = dt * (timeStageAfterImpulse(i)+1) - t_impulse_[i];
    dtau_aux_[i] = dt - dtau_[timeStageAfterImpulse(i)];
  }
  for (int i=0; i<num_lift_stages_; ++i) {
    dtau_[timeStageAfterLift(i)] = dt * (timeStageAfterLift(i)+1) - t_lift_[i];
    dtau_lift_[i] = dt - dtau_[timeStageAfterLift(i)];
  }
}

} // namespace idocp

#endif // IDOCP_PARNMPC_DISCRETIZER_HXX_ 