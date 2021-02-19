#ifndef IDOCP_OCP_DISCRETIZER_HXX_
#define IDOCP_OCP_DISCRETIZER_HXX_

#include "idocp/hybrid/ocp_discretizer.hpp"

#include <cmath>
#include <cassert>
#include <stdexcept>

namespace idocp {

inline OCPDiscretizer::OCPDiscretizer(const double T, const int N, 
                                      const int max_events, 
                                      const double sampling_period) 
  : T_(T),
    sampling_period_(sampling_period),
    N_(N),
    N_ideal_(N),
    max_events_(max_events),
    num_impulse_stages_(0),
    num_lift_stages_(0),
    contact_phase_index_from_time_stage_(N+1, 0), 
    impulse_index_after_time_stage_(N+1, -1), 
    lift_index_after_time_stage_(N+1, -1), 
    time_stage_before_impulse_(max_events, -1), 
    time_stage_before_lift_(max_events, -1),
    is_time_stage_before_impulse_(N+1, false),
    is_time_stage_before_lift_(N+1, false),
    t_(N+1, 0),
    t_impulse_(max_events, 0),
    t_lift_(max_events, 0),
    dtau_(N+1, (double)(T/N)),
    dtau_aux_(max_events, 0),
    dtau_lift_(max_events, 0) {
  try {
    if (sampling_period < 0) {
      throw std::out_of_range(
          "invalid value: sampling_period must be non-negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline OCPDiscretizer::OCPDiscretizer()
  : T_(0),
    sampling_period_(0),
    N_(0),
    max_events_(0),
    num_impulse_stages_(0),
    num_lift_stages_(0),
    contact_phase_index_from_time_stage_(), 
    impulse_index_after_time_stage_(), 
    lift_index_after_time_stage_(), 
    time_stage_before_impulse_(), 
    time_stage_before_lift_(),
    is_time_stage_before_impulse_(),
    is_time_stage_before_lift_(),
    t_(),
    t_impulse_(),
    t_lift_(),
    dtau_(),
    dtau_aux_(),
    dtau_lift_() {
}


inline OCPDiscretizer::~OCPDiscretizer() {
}


inline void OCPDiscretizer::discretizeOCP(
    const ContactSequence& contact_sequence, const double t) {
  countImpulseEvents(contact_sequence, t);
  countLiftEvents(contact_sequence, t);
  countIsTimeStageBeforeEvents(contact_sequence);
  countContactPhase(contact_sequence);
  countTime(contact_sequence, t);
}


inline int OCPDiscretizer::N() const {
  return N_;
}


inline int OCPDiscretizer::numImpulseStages() const {
  return num_impulse_stages_;
}


inline int OCPDiscretizer::numLiftStages() const {
  return num_lift_stages_;
}


inline int OCPDiscretizer::contactPhase(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N_);
  return contact_phase_index_from_time_stage_[time_stage];
}


inline int OCPDiscretizer::impulseIndexAfterTimeStage(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N_);
  return impulse_index_after_time_stage_[time_stage];
}


inline int OCPDiscretizer::liftIndexAfterTimeStage(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N_);
  return lift_index_after_time_stage_[time_stage];
}


inline int OCPDiscretizer::timeStageBeforeImpulse(
    const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < numImpulseStages());
  return time_stage_before_impulse_[impulse_index];
}


inline int OCPDiscretizer::timeStageAfterImpulse(
    const int impulse_index) const {
  return timeStageBeforeImpulse(impulse_index) + 1;
}


inline int OCPDiscretizer::timeStageBeforeLift(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < numLiftStages());
  return time_stage_before_lift_[lift_index];
}


inline int OCPDiscretizer::timeStageAfterLift(const int lift_index) const {
  return timeStageBeforeLift(lift_index) + 1;
}


inline int OCPDiscretizer::contactPhaseBeforeImpulse(
    const int impulse_index) const {
  return contactPhase(timeStageBeforeImpulse(impulse_index));
}


inline int OCPDiscretizer::contactPhaseAfterImpulse(
    const int impulse_index) const {
  return contactPhase(timeStageAfterImpulse(impulse_index));
}


inline int OCPDiscretizer::contactPhaseBeforeLift(const int lift_index) const {
  return contactPhase(timeStageBeforeLift(lift_index));
}


inline int OCPDiscretizer::contactPhaseAfterLift(const int lift_index) const {
  return contactPhase(timeStageAfterLift(lift_index));
}


inline bool OCPDiscretizer::isTimeStageBeforeImpulse(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N_);
  return is_time_stage_before_impulse_[time_stage];
}


inline bool OCPDiscretizer::isTimeStageAfterImpulse(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N_);
  return isTimeStageBeforeImpulse(time_stage-1);
}


inline bool OCPDiscretizer::isTimeStageBeforeLift(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N_);
  return is_time_stage_before_lift_[time_stage];
}


inline bool OCPDiscretizer::isTimeStageAfterLift(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N_);
  return isTimeStageBeforeLift(time_stage-1);
}


inline double OCPDiscretizer::t(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N_);
  return t_[time_stage];
}


inline double OCPDiscretizer::t_impulse(const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < numImpulseStages());
  return t_impulse_[impulse_index];
}


inline double OCPDiscretizer::t_lift(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < numLiftStages());
  return t_lift_[lift_index];
}


inline double OCPDiscretizer::dtau(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N_);
  return dtau_[time_stage];
}


inline double OCPDiscretizer::dtau_aux(const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < numImpulseStages());
  return dtau_aux_[impulse_index];
}


inline double OCPDiscretizer::dtau_lift(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < numLiftStages());
  return dtau_lift_[lift_index];
}


inline void OCPDiscretizer::countImpulseEvents(
    const ContactSequence& contact_sequence, const double t) {
  for (auto& e : impulse_index_after_time_stage_) { e = -1; }
  for (auto& e : time_stage_before_impulse_)      { e = -1; }
  num_impulse_stages_ = contact_sequence.numImpulseEvents();
  assert(num_impulse_stages_ <= max_events_);
  const double dt = T_ / N_;
  for (int impulse_index=0; impulse_index<num_impulse_stages_; ++impulse_index) {
    t_impulse_[impulse_index] = contact_sequence.impulseTime(impulse_index) - t;
    time_stage_before_impulse_[impulse_index] 
        = std::floor(t_impulse_[impulse_index]/dt);
    impulse_index_after_time_stage_[time_stage_before_impulse_[impulse_index]] 
        = impulse_index;
  }
}


inline void OCPDiscretizer::countLiftEvents(
    const ContactSequence& contact_sequence, const double t) {
  for (auto& e : lift_index_after_time_stage_) { e = -1; }
  for (auto& e : time_stage_before_lift_)      { e = -1; }
  num_lift_stages_ = contact_sequence.numLiftEvents();
  assert(num_lift_stages_ <=  max_events_);
  const double dt = T_ / N_;
  for (int lift_index=0; lift_index<num_lift_stages_; ++lift_index) {
    t_lift_[lift_index] = contact_sequence.liftTime(lift_index) - t;
    time_stage_before_lift_[lift_index] = std::floor(t_lift_[lift_index]/dt);
    lift_index_after_time_stage_[time_stage_before_lift_[lift_index]] 
        = lift_index;
  }
}


inline void OCPDiscretizer::countIsTimeStageBeforeEvents(
    const ContactSequence& contact_sequence) {
  int impulse_index = 0;
  int lift_index = 0;
  for (int i=0; i<N_; ++i) {
    if (impulse_index < numImpulseStages()) {
      if (i == timeStageBeforeImpulse(impulse_index)) {
        is_time_stage_before_impulse_[i] = true;
        ++impulse_index;
      }
      else {
        is_time_stage_before_impulse_[i] = false;
      }
    }
    else {
      is_time_stage_before_impulse_[i] = false;
    }
    if (lift_index < numLiftStages()) {
      if (i == timeStageBeforeLift(lift_index)) {
        is_time_stage_before_lift_[i] = true;
        ++lift_index;
      }
      else {
        is_time_stage_before_lift_[i] = false;
      }
    }
    else {
      is_time_stage_before_lift_[i] = false;
    }
  }
}


inline void OCPDiscretizer::countContactPhase(
    const ContactSequence& contact_sequence) {
  for (auto & e : contact_phase_index_from_time_stage_) { e = 0; }
  int num_events = 0;
  for (int i=0; i<=N_; ++i) {
    contact_phase_index_from_time_stage_[i] = num_events;
    if (isTimeStageBeforeImpulse(i) || isTimeStageBeforeLift(i)) {
      ++num_events; 
    }
  }
}


inline void OCPDiscretizer::countTime(const ContactSequence& contact_sequence, 
                                      const double t) {
  const double dt = T_ / N_;
  for (int i=0; i<=N_; ++i) {
    t_[i] = t + i * dt;
    dtau_[i] = dt;
  }
  for (int i=0; i<num_impulse_stages_; ++i) {
    dtau_[timeStageBeforeImpulse(i)] 
        = t_impulse_[i] - dt * timeStageBeforeImpulse(i);
    dtau_aux_[i] = dt - dtau_[timeStageBeforeImpulse(i)];
  }
  for (int i=0; i<num_lift_stages_; ++i) {
    dtau_[timeStageBeforeLift(i)] = t_lift_[i] - dt * timeStageBeforeLift(i);
    dtau_lift_[i] = dt - dtau_[timeStageBeforeLift(i)];
  }
}

} // namespace idocp

#endif // IDOCP_OCP_DISCRETIZER_HXX_ 