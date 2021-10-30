#ifndef ROBOTOC_HYBRID_OCP_DISCRETIZATION_HXX_
#define ROBOTOC_HYBRID_OCP_DISCRETIZATION_HXX_ 

#include "robotoc/hybrid/hybrid_ocp_discretization.hpp"

#include <cassert>

namespace robotoc {

inline HybridOCPDiscretization::HybridOCPDiscretization(const double T, 
                                                        const int N, 
                                                        const int max_events) 
  : T_(T),
    dt_ideal_(T/N), 
    max_dt_(dt_ideal_-k_min_dt),
    N_(N),
    N_ideal_(N),
    N_impulse_(0),
    N_lift_(0),
    max_events_(max_events),
    contact_phase_index_from_time_stage_(N+1, 0), 
    impulse_index_after_time_stage_(N+1, -1), 
    lift_index_after_time_stage_(N+1, -1), 
    time_stage_before_impulse_(max_events+1, -1), 
    time_stage_before_lift_(max_events+1, -1),
    is_time_stage_before_impulse_(N+1, false),
    is_time_stage_before_lift_(N+1, false),
    t_(N+1, 0),
    t_impulse_(max_events+1, 0),
    t_lift_(max_events+1, 0),
    dt_(N+1, static_cast<double>(T/N)),
    dt_aux_(max_events+1, 0),
    dt_lift_(max_events+1, 0),
    event_types_(2*max_events+1, DiscreteEventType::None),
    sto_impulse_(max_events), 
    sto_lift_(max_events) {
}


inline HybridOCPDiscretization::HybridOCPDiscretization()
  : T_(0),
    dt_ideal_(0), 
    max_dt_(0),
    N_(0),
    N_ideal_(0),
    N_impulse_(0),
    N_lift_(0),
    max_events_(0),
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
    dt_(),
    dt_aux_(),
    dt_lift_(),
    event_types_(),
    sto_impulse_(), 
    sto_lift_() {
}


inline HybridOCPDiscretization::~HybridOCPDiscretization() {
}


inline void HybridOCPDiscretization::discretize(
    const std::shared_ptr<ContactSequence>& contact_sequence, const double t) {
  countDiscreteEvents(contact_sequence, t);
  countTimeSteps(t);
  countTimeStages();
  countContactPhase();
  assert(isFormulationTractable());
  assert(isSwitchingTimeConsistent());
}


inline int HybridOCPDiscretization::N() const {
  return N_;
}


inline int HybridOCPDiscretization::N_impulse() const {
  return N_impulse_;
}


inline int HybridOCPDiscretization::N_lift() const {
  return N_lift_;
}


inline int HybridOCPDiscretization::N_ideal() const {
  return N_ideal_;
}


inline int HybridOCPDiscretization::contactPhase(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N());
  return contact_phase_index_from_time_stage_[time_stage];
}


inline int HybridOCPDiscretization::contactPhaseAfterImpulse(
    const int impulse_index) const {
  return contactPhase(timeStageAfterImpulse(impulse_index));
}


inline int HybridOCPDiscretization::contactPhaseAfterLift(
    const int lift_index) const {
  return contactPhase(timeStageAfterLift(lift_index));
}


inline int HybridOCPDiscretization::impulseIndexAfterTimeStage(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N());
  return impulse_index_after_time_stage_[time_stage];
}


inline int HybridOCPDiscretization::liftIndexAfterTimeStage(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N());
  return lift_index_after_time_stage_[time_stage];
}


inline int HybridOCPDiscretization::timeStageBeforeImpulse(
    const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < N_impulse());
  return time_stage_before_impulse_[impulse_index];
}


inline int HybridOCPDiscretization::timeStageAfterImpulse(
    const int impulse_index) const {
  return timeStageBeforeImpulse(impulse_index) + 1;
}


inline int HybridOCPDiscretization::timeStageBeforeLift(
    const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < N_lift());
  return time_stage_before_lift_[lift_index];
}


inline int HybridOCPDiscretization::timeStageAfterLift(
    const int lift_index) const {
  return timeStageBeforeLift(lift_index) + 1;
}


inline bool HybridOCPDiscretization::isTimeStageBeforeImpulse(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N());
  return is_time_stage_before_impulse_[time_stage];
}


inline bool HybridOCPDiscretization::isTimeStageAfterImpulse(
    const int time_stage) const {
  assert(time_stage > 0);
  assert(time_stage <= N());
  return isTimeStageBeforeImpulse(time_stage-1);
}


inline bool HybridOCPDiscretization::isTimeStageBeforeLift(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N());
  return is_time_stage_before_lift_[time_stage];
}


inline bool HybridOCPDiscretization::isTimeStageAfterLift(
    const int time_stage) const {
  assert(time_stage > 0);
  assert(time_stage <= N());
  return isTimeStageBeforeLift(time_stage-1);
}


inline double HybridOCPDiscretization::t(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N());
  return t_[time_stage];
}


inline double HybridOCPDiscretization::t_impulse(
    const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < N_impulse());
  return t_impulse_[impulse_index];
}


inline double HybridOCPDiscretization::t_lift(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < N_lift());
  return t_lift_[lift_index];
}


inline double HybridOCPDiscretization::dt(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N());
  return dt_[time_stage];
}


inline double HybridOCPDiscretization::dt_aux(const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < N_impulse());
  return dt_aux_[impulse_index];
}


inline double HybridOCPDiscretization::dt_lift(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < N_lift());
  return dt_lift_[lift_index];
}


inline double HybridOCPDiscretization::dt_ideal() const {
  return dt_ideal_;
}


inline bool HybridOCPDiscretization::isSTOEnabledImpulse(
    const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < N_impulse());
  return sto_impulse_[impulse_index];
}


inline bool HybridOCPDiscretization::isSTOEnabledLift(
    const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < N_lift());
  return sto_lift_[lift_index];
}


inline int HybridOCPDiscretization::eventIndexImpulse(
    const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < N_impulse());
  return (contactPhaseAfterImpulse(impulse_index)-1);
}


inline int HybridOCPDiscretization::eventIndexLift(
    const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < N_lift());
  return (contactPhaseAfterLift(lift_index)-1);
}


inline DiscreteEventType HybridOCPDiscretization::eventType(
    const int event_index) const {
  assert(event_index < N_impulse()+N_lift());
  return event_types_[event_index];
}


inline bool HybridOCPDiscretization::isFormulationTractable() const {
  for (int i=0; i<N(); ++i) {
    if (isTimeStageBeforeImpulse(i) && isTimeStageBeforeLift(i)) {
      return false;
    }
  }
  for (int i=0; i<N()-1; ++i) {
    if (isTimeStageBeforeImpulse(i) && isTimeStageBeforeImpulse(i+1)) {
      return false;
    }
    if (isTimeStageBeforeLift(i) && isTimeStageBeforeImpulse(i+1)) {
      return false;
    }
  }
  return true;
}


inline bool HybridOCPDiscretization::isSwitchingTimeConsistent() const {
  for (int i=0; i<N_impulse(); ++i) {
    if (t_impulse(i) < t(0)+k_min_dt || t_impulse(i) >= t(N())-k_min_dt) {
      return false;
    }
  }
  for (int i=0; i<N_lift(); ++i) {
    if (t_lift(i) < t(0)+k_min_dt || t_lift(i) > t(N())-k_min_dt) {
      return false;
    }
  }
  return true;
}


inline void HybridOCPDiscretization::countDiscreteEvents(
    const std::shared_ptr<ContactSequence>& contact_sequence, const double t) {
  const int max_num_impulse_events = contact_sequence->numImpulseEvents();
  assert(max_num_impulse_events <= max_events_);
  N_impulse_ = 0;
  for (int impulse_index=0; impulse_index<max_num_impulse_events; ++impulse_index) {
    const double t_impulse = contact_sequence->impulseTime(impulse_index);
    if (t_impulse >= t+T_-k_min_dt) {
      break;
    }
    t_impulse_[impulse_index] = t_impulse;
    time_stage_before_impulse_[impulse_index] = std::floor((t_impulse-t)/dt_ideal_);
    sto_impulse_[impulse_index] = contact_sequence->isSTOEnabledImpulse(impulse_index);
    ++N_impulse_;
  }
  const int max_num_lift_events = contact_sequence->numLiftEvents();
  assert(max_num_lift_events <= max_events_);
  N_lift_ = 0;
  for (int lift_index=0; lift_index<max_num_lift_events; ++lift_index) {
    const double t_lift = contact_sequence->liftTime(lift_index);
    if (t_lift >= t+T_-k_min_dt) {
      break;
    }
    t_lift_[lift_index] = t_lift;
    time_stage_before_lift_[lift_index] = std::floor((t_lift-t)/dt_ideal_);
    sto_lift_[lift_index] = contact_sequence->isSTOEnabledLift(lift_index);
    ++N_lift_;
  }
  N_ = N_ideal_;
  for (int i=0; i<N_impulse_+N_lift_; ++i) {
    event_types_[i] = contact_sequence->eventType(i);
  } 
  assert(isFormulationTractable());
}


inline void HybridOCPDiscretization::countTimeSteps(const double t) {
  int impulse_index = 0;
  int lift_index = 0;
  int num_events_on_grid = 0;
  for (int i=0; i<N_ideal_; ++i) {
    const int stage = i - num_events_on_grid;
    if (i == time_stage_before_impulse_[impulse_index]) {
      dt_[stage] = t_impulse_[impulse_index] - i * dt_ideal_ - t;
      assert(dt_[stage] >= -k_min_dt);
      assert(dt_[stage] <= dt_ideal_+k_min_dt);
      if (dt_[stage] <= k_min_dt) {
        time_stage_before_impulse_[impulse_index] = stage - 1;
        dt_aux_[impulse_index] = dt_ideal_;
        t_[stage] = t + (i-1) * dt_ideal_;
        ++num_events_on_grid;
        ++impulse_index;
      }
      else if (dt_[stage] >= max_dt_) {
        time_stage_before_impulse_[impulse_index] = i + 1;
        t_[stage] = t + i * dt_ideal_;
      }
      else {
        time_stage_before_impulse_[impulse_index] = stage;
        dt_aux_[impulse_index] = dt_ideal_ - dt_[stage];
        t_[stage] = t + i * dt_ideal_;
        ++impulse_index;
      }
    }
    else if (i == time_stage_before_lift_[lift_index]) {
      dt_[stage] = t_lift_[lift_index] - i * dt_ideal_ - t;
      assert(dt_[stage] >= -k_min_dt);
      assert(dt_[stage] <= dt_ideal_+k_min_dt);
      if (dt_[stage] <= k_min_dt) {
        time_stage_before_lift_[lift_index] = stage - 1;
        dt_lift_[lift_index] = dt_ideal_;
        t_[stage] = t + (i-1) * dt_ideal_;
        ++num_events_on_grid;
        ++lift_index;
      }
      else if (dt_[stage] >= max_dt_) {
        time_stage_before_lift_[lift_index] = i + 1;
        t_[stage] = t + i * dt_ideal_;
      }
      else {
        time_stage_before_lift_[lift_index] = stage;
        dt_lift_[lift_index] = dt_ideal_ - dt_[stage];
        t_[stage] = t + i * dt_ideal_;
        ++lift_index;
      }
    }
    else {
      dt_[stage] = dt_ideal_;
      t_[stage] = t + i * dt_ideal_;
    }
  }
  N_ = N_ideal_ - num_events_on_grid;
  t_[N_] = t + T_;
}


inline void HybridOCPDiscretization::countTimeStages() {
  int impulse_index = 0;
  int lift_index = 0;
  for (int i=0; i<N(); ++i) {
    if (impulse_index < N_impulse()) {
      if (i == timeStageBeforeImpulse(impulse_index)) {
        is_time_stage_before_impulse_[i] = true;
        impulse_index_after_time_stage_[i] = impulse_index; 
        ++impulse_index;
      }
      else {
        is_time_stage_before_impulse_[i] = false;
        impulse_index_after_time_stage_[i] = -1; 
      }
    }
    else {
      is_time_stage_before_impulse_[i] = false;
      impulse_index_after_time_stage_[i] = -1; 
    }
    if (lift_index < N_lift()) {
      if (i == timeStageBeforeLift(lift_index)) {
        is_time_stage_before_lift_[i] = true;
        lift_index_after_time_stage_[i] = lift_index;
        ++lift_index;
      }
      else {
        is_time_stage_before_lift_[i] = false;
        lift_index_after_time_stage_[i] = -1;
      }
    }
    else {
      is_time_stage_before_lift_[i] = false;
      lift_index_after_time_stage_[i] = -1;
    }
  }
  is_time_stage_before_impulse_[N()] = false;
  is_time_stage_before_lift_[N()] = false;
}


inline void HybridOCPDiscretization::countContactPhase() {
  int num_events = 0;
  for (int i=0; i<N(); ++i) {
    contact_phase_index_from_time_stage_[i] = num_events;
    if (isTimeStageBeforeImpulse(i) || isTimeStageBeforeLift(i)) {
      ++num_events; 
    }
  }
  contact_phase_index_from_time_stage_[N()] = num_events;
}

} // namespace robotoc

#endif // ROBOTOC_HYBRID_OCP_DISCRETIZATION_HXX_ 