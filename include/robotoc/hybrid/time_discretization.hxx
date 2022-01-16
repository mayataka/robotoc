#ifndef ROBOTOC_TIME_DISCRETIZATION_HXX_
#define ROBOTOC_TIME_DISCRETIZATION_HXX_ 

#include "robotoc/hybrid/time_discretization.hpp"

#include <cassert>
#include <numeric>


namespace robotoc {

inline TimeDiscretization::TimeDiscretization(const double T, const int N, 
                                              const int max_num_each_discrete_events) 
  : T_(T),
    dt_ideal_(T/N), 
    max_dt_(dt_ideal_-k_min_dt),
    N_(N),
    N_ideal_(N),
    N_impulse_(0),
    N_lift_(0),
    max_num_each_discrete_events_(max_num_each_discrete_events),
    N_phase_(2*max_num_each_discrete_events+1, 1),
    contact_phase_from_time_stage_(N+1, 0), 
    impulse_index_after_time_stage_(N+1, -1), 
    lift_index_after_time_stage_(N+1, -1), 
    time_stage_before_impulse_(max_num_each_discrete_events+1, -1), 
    time_stage_before_lift_(max_num_each_discrete_events+1, -1),
    is_time_stage_before_impulse_(N+1, false),
    is_time_stage_before_lift_(N+1, false),
    grid_(N+1, GridInfo()), 
    grid_impulse_(max_num_each_discrete_events+1, GridInfo()), 
    grid_lift_(max_num_each_discrete_events+1, GridInfo()),
    event_types_(2*max_num_each_discrete_events+1, DiscreteEventType::None),
    sto_impulse_(max_num_each_discrete_events), 
    sto_lift_(max_num_each_discrete_events),
    sto_event_(2*max_num_each_discrete_events+1),
    discretization_method_(DiscretizationMethod::GridBased) {
}


inline TimeDiscretization::TimeDiscretization()
  : T_(0),
    dt_ideal_(0), 
    max_dt_(0),
    N_(0),
    N_ideal_(0),
    N_impulse_(0),
    N_lift_(0),
    max_num_each_discrete_events_(0),
    N_phase_(),
    contact_phase_from_time_stage_(), 
    impulse_index_after_time_stage_(), 
    lift_index_after_time_stage_(), 
    time_stage_before_impulse_(), 
    time_stage_before_lift_(),
    is_time_stage_before_impulse_(),
    is_time_stage_before_lift_(),
    grid_(), 
    grid_impulse_(), 
    grid_lift_(),
    event_types_(),
    sto_impulse_(), 
    sto_lift_(),
    sto_event_(),
    discretization_method_(DiscretizationMethod::GridBased) {
}


inline TimeDiscretization::~TimeDiscretization() {
}


inline void TimeDiscretization::setDiscretizationMethod(
    const DiscretizationMethod discretization_method) {
  discretization_method_ = discretization_method;
}


inline void TimeDiscretization::discretize(
    const std::shared_ptr<ContactSequence>& contact_sequence, const double t) {
  if (discretization_method_ == DiscretizationMethod::GridBased) {
    countDiscreteEvents(contact_sequence, t, true);
    countTimeStepsGridBased(t);
  }
  else {
    countDiscreteEvents(contact_sequence, t, false);
    countTimeStepsPhaseBased(t);
  }
  countTimeStages();
  countContactPhase();
  countSTOEvents();
  assert(isFormulationTractable());
  assert(isSwitchingTimeConsistent());
}


inline void TimeDiscretization::meshRefinement(
    const std::shared_ptr<ContactSequence>& contact_sequence, const double t) {
  if (discretization_method_ == DiscretizationMethod::PhaseBased) {
    countDiscreteEvents(contact_sequence, t, true);
    countTimeStepsPhaseBased(t);
    countTimeStages();
    countContactPhase();
    countSTOEvents();
    assert(isFormulationTractable());
    assert(isSwitchingTimeConsistent());
  }
}


inline int TimeDiscretization::N() const {
  return N_;
}


inline int TimeDiscretization::N_impulse() const {
  return N_impulse_;
}


inline int TimeDiscretization::N_lift() const {
  return N_lift_;
}


inline int TimeDiscretization::N_ideal() const {
  return N_ideal_;
}


inline int TimeDiscretization::N_phase(const int phase) const {
  assert(phase >= 0);
  assert(phase <= N_impulse()+N_lift());
  return N_phase_[phase];
}


inline int TimeDiscretization::numContactPhases() const {
  return (numDiscreteEvents()+1);
}


inline int TimeDiscretization::numDiscreteEvents() const {
  return (N_impulse_+N_lift_);
}


inline int TimeDiscretization::contactPhase(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N());
  return grid_[time_stage].contact_phase;
}


inline int TimeDiscretization::contactPhaseAfterImpulse(
    const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < N_impulse());
  return grid_impulse_[impulse_index].contact_phase;
}


inline int TimeDiscretization::contactPhaseAfterLift(
    const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < N_lift());
  return grid_lift_[lift_index].contact_phase;
}


inline int TimeDiscretization::impulseIndexAfterTimeStage(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N());
  return impulse_index_after_time_stage_[time_stage];
}


inline int TimeDiscretization::liftIndexAfterTimeStage(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N());
  return lift_index_after_time_stage_[time_stage];
}


inline int TimeDiscretization::timeStageBeforeImpulse(
    const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < N_impulse());
  return time_stage_before_impulse_[impulse_index];
}


inline int TimeDiscretization::timeStageAfterImpulse(
    const int impulse_index) const {
  return timeStageBeforeImpulse(impulse_index) + 1;
}


inline int TimeDiscretization::timeStageBeforeLift(
    const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < N_lift());
  return time_stage_before_lift_[lift_index];
}


inline int TimeDiscretization::timeStageAfterLift(
    const int lift_index) const {
  return timeStageBeforeLift(lift_index) + 1;
}


inline bool TimeDiscretization::isTimeStageBeforeImpulse(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N());
  return is_time_stage_before_impulse_[time_stage];
}


inline bool TimeDiscretization::isTimeStageAfterImpulse(
    const int time_stage) const {
  assert(time_stage > 0);
  assert(time_stage <= N());
  return isTimeStageBeforeImpulse(time_stage-1);
}


inline bool TimeDiscretization::isTimeStageBeforeLift(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N());
  return is_time_stage_before_lift_[time_stage];
}


inline bool TimeDiscretization::isTimeStageAfterLift(
    const int time_stage) const {
  assert(time_stage > 0);
  assert(time_stage <= N());
  return isTimeStageBeforeLift(time_stage-1);
}


inline double TimeDiscretization::t(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N());
  return grid_[time_stage].t;
}


inline double TimeDiscretization::t_impulse(const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < N_impulse());
  return grid_impulse_[impulse_index].t;
}


inline double TimeDiscretization::t_lift(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < N_lift());
  return grid_lift_[lift_index].t;
}


inline double TimeDiscretization::dt(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N());
  return grid_[time_stage].dt;
}


inline double TimeDiscretization::dt_aux(const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < N_impulse());
  return grid_impulse_[impulse_index].dt;
}


inline double TimeDiscretization::dt_lift(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < N_lift());
  return grid_lift_[lift_index].dt;
}


inline double TimeDiscretization::dt_max() const {
  std::vector<double> dt_phase;
  dt_phase.push_back(dt(0));
  for (int impulse_index=0; impulse_index<N_impulse(); ++impulse_index) {
    dt_phase.push_back(dt_aux(impulse_index));
  }
  for (int lift_index=0; lift_index<N_lift(); ++lift_index) {
    dt_phase.push_back(dt_lift(lift_index));
  }
  return *std::max_element(dt_phase.begin(), dt_phase.end());
}


inline double TimeDiscretization::dt_ideal() const {
  return dt_ideal_;
}


inline const GridInfo& TimeDiscretization::gridInfo(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N_);
  return grid_[time_stage];
}


inline const GridInfo& TimeDiscretization::gridInfoImpulse(
    const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < N_impulse());
  return grid_impulse_[impulse_index];
}


inline const GridInfo& TimeDiscretization::gridInfoAux(
    const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < N_impulse());
  return grid_impulse_[impulse_index];
}


inline const GridInfo& TimeDiscretization::gridInfoLift(
    const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < N_impulse());
  return grid_lift_[lift_index];
}


inline int TimeDiscretization::timeStageInPhase(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N_);
  return grid_[time_stage].grid_count_in_phase;
}


inline bool TimeDiscretization::isSTOEnabledEvent(
    const int event_index) const {
  assert(event_index >= 0);
  assert(event_index < N_impulse()+N_lift());
  return sto_event_[event_index];
}


inline bool TimeDiscretization::isSTOEnabledPhase(const int phase) const {
  assert(phase >= 0);
  assert(phase < numContactPhases());
  if (phase == 0) {
    return sto_event_[0];
  }
  else if (phase == numContactPhases()-1) {
    return sto_event_[numDiscreteEvents()-1];
  }
  else {
    return (sto_event_[phase-1] || sto_event_[phase]);
  }
}


inline bool TimeDiscretization::isSTOEnabledNextPhase(const int phase) const {
  if (phase+1 == numContactPhases()) return false;
  else return isSTOEnabledPhase(phase+1);
}


inline bool TimeDiscretization::isSTOEnabledImpulse(
    const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < N_impulse());
  return sto_impulse_[impulse_index];
}


inline bool TimeDiscretization::isSTOEnabledLift(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < N_lift());
  return sto_lift_[lift_index];
}


inline int TimeDiscretization::eventIndexImpulse(
    const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < N_impulse());
  return (contactPhaseAfterImpulse(impulse_index)-1);
}


inline int TimeDiscretization::eventIndexLift(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < N_lift());
  return (contactPhaseAfterLift(lift_index)-1);
}


inline DiscreteEventType TimeDiscretization::eventType(
    const int event_index) const {
  assert(event_index < N_impulse()+N_lift());
  return event_types_[event_index];
}


inline DiscretizationMethod TimeDiscretization::discretizationMethod() const {
  return discretization_method_;
}


inline int TimeDiscretization::maxNumEachDiscreteEvents() const {
  return max_num_each_discrete_events_;
}


inline std::vector<double> TimeDiscretization::timeSteps() const {
  std::vector<double> time_steps;
  for (int i=0; i<N(); ++i) {
    time_steps.push_back(dt(i));
    if (isTimeStageBeforeImpulse(i)) {
      time_steps.push_back(dt_aux(impulseIndexAfterTimeStage(i)));
    }
    else if (isTimeStageBeforeLift(i)) {
      time_steps.push_back(dt_lift(liftIndexAfterTimeStage(i)));
    }
  }
  return time_steps;
}


inline std::vector<double> TimeDiscretization::timePoints() const {
  std::vector<double> time_points;
  for (int i=0; i<N(); ++i) {
    time_points.push_back(t(i));
    if (isTimeStageBeforeImpulse(i)) {
      time_points.push_back(t_impulse(impulseIndexAfterTimeStage(i)));
    }
    else if (isTimeStageBeforeLift(i)) {
      time_points.push_back(t_lift(liftIndexAfterTimeStage(i)));
    }
  }
  time_points.push_back(t(N()));
  return time_points;
}


inline bool TimeDiscretization::isFormulationTractable() const {
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


inline bool TimeDiscretization::isSwitchingTimeConsistent() const {
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


inline void TimeDiscretization::countDiscreteEvents(
    const std::shared_ptr<ContactSequence>& contact_sequence, const double t,
    const bool refine_grids) {
  const int max_num_impulse_events = contact_sequence->numImpulseEvents();
  assert(max_num_impulse_events <= max_num_each_discrete_events_);
  const bool is_phase_based 
    = (discretization_method_ == DiscretizationMethod::PhaseBased);
  N_impulse_ = 0;
  for (int impulse_index=0; impulse_index<max_num_impulse_events; ++impulse_index) {
    const double t_impulse = contact_sequence->impulseTime(impulse_index);
    if (t_impulse >= t+T_-k_min_dt) {
      break;
    }
    grid_impulse_[impulse_index].t = t_impulse;
    if (refine_grids) {
      time_stage_before_impulse_[impulse_index] = std::floor((t_impulse-t)/dt_ideal_);
    }
    sto_impulse_[impulse_index] 
        = (is_phase_based && contact_sequence->isSTOEnabledImpulse(impulse_index));
    ++N_impulse_;
  }
  const int max_num_lift_events = contact_sequence->numLiftEvents();
  assert(max_num_lift_events <= max_num_each_discrete_events_);
  N_lift_ = 0;
  for (int lift_index=0; lift_index<max_num_lift_events; ++lift_index) {
    const double t_lift = contact_sequence->liftTime(lift_index);
    if (t_lift >= t+T_-k_min_dt) {
      break;
    }
    grid_lift_[lift_index].t = t_lift;
    if (refine_grids) {
      time_stage_before_lift_[lift_index] = std::floor((t_lift-t)/dt_ideal_);
    }
    sto_lift_[lift_index] 
        = (is_phase_based && contact_sequence->isSTOEnabledLift(lift_index));
    ++N_lift_;
  }
  N_ = N_ideal_;
  for (int i=0; i<N_impulse_+N_lift_; ++i) {
    event_types_[i] = contact_sequence->eventType(i);
  } 
  assert(isFormulationTractable());
}


inline void TimeDiscretization::countTimeStepsGridBased(const double t) {
  int impulse_index = 0;
  int lift_index = 0;
  int num_events_on_grid = 0;
  int time_stage_before_prev_event = 0;
  int phase = 0;
  for (int i=0; i<N_ideal_; ++i) {
    const int stage = i - num_events_on_grid;
    if (i == time_stage_before_impulse_[impulse_index]) {
      grid_[stage].dt = grid_impulse_[impulse_index].t - i * dt_ideal_ - t;
      assert(grid_[stage].dt >= -k_min_dt);
      assert(grid_[stage].dt <= dt_ideal_+k_min_dt);
      if (grid_[stage].dt <= k_min_dt) {
        time_stage_before_impulse_[impulse_index] = stage - 1;
        grid_impulse_[impulse_index].dt = dt_ideal_;
        grid_[stage].t = t + (i-1) * dt_ideal_;
        ++num_events_on_grid;
        ++impulse_index;
      }
      else if (grid_[stage].dt >= max_dt_) {
        time_stage_before_impulse_[impulse_index] = i + 1;
        grid_[stage].t = t + i * dt_ideal_;
      }
      else {
        time_stage_before_impulse_[impulse_index] = stage;
        grid_impulse_[impulse_index].dt = dt_ideal_ - grid_[stage].dt;
        grid_[stage].t = t + i * dt_ideal_;
        ++impulse_index;
      }
      N_phase_[phase] = time_stage_before_impulse_[impulse_index-1] 
                          - time_stage_before_prev_event + 1;
      time_stage_before_prev_event = time_stage_before_impulse_[impulse_index-1];
      ++phase;
    }
    else if (i == time_stage_before_lift_[lift_index]) {
      grid_[stage].dt = grid_lift_[lift_index].t - i * dt_ideal_ - t;
      assert(grid_[stage].dt >= -k_min_dt);
      assert(grid_[stage].dt <= dt_ideal_+k_min_dt);
      if (grid_[stage].dt <= k_min_dt) {
        time_stage_before_lift_[lift_index] = stage - 1;
        grid_lift_[lift_index].dt = dt_ideal_;
        grid_[stage].t = t + (i-1) * dt_ideal_;
        ++num_events_on_grid;
        ++lift_index;
      }
      else if (grid_[stage].dt >= max_dt_) {
        time_stage_before_lift_[lift_index] = i + 1;
        grid_[stage].t = t + i * dt_ideal_;
      }
      else {
        time_stage_before_lift_[lift_index] = stage;
        grid_lift_[lift_index].dt = dt_ideal_ - grid_[stage].dt;
        grid_[stage].t = t + i * dt_ideal_;
        ++lift_index;
      }
      N_phase_[phase] = time_stage_before_lift_[lift_index-1] 
                          - time_stage_before_prev_event + 1;
      time_stage_before_prev_event = time_stage_before_lift_[lift_index-1];
      ++phase;
    }
    else {
      grid_[stage].dt = dt_ideal_;
      grid_[stage].t = t + i * dt_ideal_;
    }
    grid_[stage].grid_count_in_phase = stage - time_stage_before_prev_event;
  }
  N_ = N_ideal_ - num_events_on_grid;
  grid_[N_].t = t + T_;
  N_phase_[phase] = N_ - time_stage_before_prev_event;
  grid_[N_].grid_count_in_phase = N_ - time_stage_before_prev_event;
}


inline void TimeDiscretization::countTimeStepsPhaseBased(const double t) {
  int next_impulse_index = 0;
  int next_lift_index = 0;
  int time_stage_before_prev_event = 0;
  double t_prev_event = t;
  grid_[0].t = t;
  for (int phase=0; phase<N_impulse()+N_lift(); ++phase) {
    const int next_event_index = phase;
    const auto next_event_type = eventType(next_event_index);
    assert(next_event_type != DiscreteEventType::None);
    if (next_event_type == DiscreteEventType::Impulse) {
      const int time_stage_before_next_event 
          = time_stage_before_impulse_[next_impulse_index];
      const int num_phase_grids = time_stage_before_next_event
                                  - time_stage_before_prev_event + 1;
      N_phase_[phase] = num_phase_grids;
      const double dt_phase 
          = (t_impulse(next_impulse_index)-t_prev_event) / num_phase_grids;
      for (int stage=time_stage_before_prev_event+1; 
            stage<=time_stage_before_next_event; ++stage) {
        grid_[stage].dt = dt_phase;
      }
      grid_[time_stage_before_prev_event+1].t = t_prev_event + dt_phase;
      for (int stage=time_stage_before_prev_event+2; 
            stage<=time_stage_before_next_event; ++stage) {
        grid_[stage].t = grid_[stage-1].t + dt_phase;
      }
      for (int stage=time_stage_before_prev_event+1; 
            stage<=time_stage_before_next_event; ++stage) {
        grid_[stage].grid_count_in_phase = stage - time_stage_before_prev_event;
      }
      time_stage_before_prev_event = time_stage_before_next_event;
      t_prev_event = t_impulse(next_impulse_index);
      ++next_impulse_index;
    }
    else {
      const int time_stage_before_next_event 
          = time_stage_before_lift_[next_lift_index];
      const int num_phase_grids = time_stage_before_next_event
                                  - time_stage_before_prev_event + 1;
      N_phase_[phase] = num_phase_grids;
      const double dt_phase
          = (t_lift(next_lift_index)-t_prev_event) / num_phase_grids;
      for (int stage=time_stage_before_prev_event+1; 
            stage<=time_stage_before_next_event; ++stage) {
        grid_[stage].dt = dt_phase;
      }
      grid_[time_stage_before_prev_event+1].t = t_prev_event + dt_phase;
      for (int stage=time_stage_before_prev_event+2; 
            stage<=time_stage_before_next_event; ++stage) {
        grid_[stage].t = grid_[stage-1].t + dt_phase;
      }
      for (int stage=time_stage_before_prev_event+1; 
            stage<=time_stage_before_next_event; ++stage) {
        grid_[stage].grid_count_in_phase = stage - time_stage_before_prev_event;
      }
      time_stage_before_prev_event = time_stage_before_next_event;
      t_prev_event = t_lift(next_lift_index);
      ++next_lift_index;
    }
  }
  const int last_phase = N_impulse() + N_lift();
  const int num_phase_grids = N_ - time_stage_before_prev_event;
  N_phase_[last_phase] = num_phase_grids;
  const double dt_phase = (t+T_-t_prev_event) / num_phase_grids;
  for (int stage=time_stage_before_prev_event+1; stage<N_; ++stage) {
    grid_[stage].dt = dt_phase;
  }
  grid_[time_stage_before_prev_event+1].t = t_prev_event + dt_phase;
  for (int stage=time_stage_before_prev_event+2; stage<=N_; ++stage) {
    grid_[stage].t = grid_[stage-1].t + dt_phase;
  }
  for (int stage=time_stage_before_prev_event+1; stage<=N_; ++stage) {
    grid_[stage].grid_count_in_phase = stage - time_stage_before_prev_event;
  }
  for (int impulse_index=0; impulse_index<N_impulse(); ++impulse_index) {
    grid_impulse_[impulse_index].dt 
        = grid_[time_stage_before_impulse_[impulse_index]+1].dt;
  }
  for (int lift_index=0; lift_index<N_lift(); ++lift_index) {
    grid_lift_[lift_index].dt 
        = grid_[time_stage_before_lift_[lift_index]+1].dt;
  }
  grid_[0].dt = grid_[1].dt;
  if (N_impulse() > 0) {
    if (time_stage_before_impulse_[0] == 0) {
      grid_[0].dt = grid_impulse_[0].t - t;
    }
  }
  if (N_lift() > 0) {
    if (time_stage_before_lift_[0] == 0) {
      grid_[0].dt = grid_lift_[0].t - t;
    }
  } 
}


inline void TimeDiscretization::countTimeStages() {
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


inline void TimeDiscretization::countContactPhase() {
  int num_events = 0;
  for (int i=0; i<N(); ++i) {
    grid_[i].contact_phase = num_events;
    if (isTimeStageBeforeImpulse(i) || isTimeStageBeforeLift(i)) {
      ++num_events; 
    }
  }
  grid_[N()].contact_phase = num_events;
  for (int impulse_index=0; impulse_index<N_impulse(); ++impulse_index) {
    grid_impulse_[impulse_index].contact_phase = contactPhase(timeStageAfterImpulse(impulse_index));
  }
  for (int lift_index=0; lift_index<N_lift(); ++lift_index) {
    grid_lift_[lift_index].contact_phase = contactPhase(timeStageAfterLift(lift_index));
  }
}


inline void TimeDiscretization::countSTOEvents() {
  const bool is_phase_based 
      = (discretization_method_ == DiscretizationMethod::PhaseBased);
  int event_index = 0;
  for (int i=0; i<N(); ++i) {
    if (isTimeStageBeforeImpulse(i)) {
      sto_event_[event_index] 
          = (is_phase_based && isSTOEnabledImpulse(impulseIndexAfterTimeStage(i)));
      ++event_index;
    }
    else if (isTimeStageBeforeLift(i)) {
      sto_event_[event_index] 
          = (is_phase_based && isSTOEnabledLift(liftIndexAfterTimeStage(i)));
      ++event_index;
    }
  }
  for (int i=N_lift()+N_impulse(); i<sto_event_.size(); ++i) {
    sto_event_[i] = false;
  }
}

} // namespace robotoc

#endif // ROBOTOC_TIME_DISCRETIZATION_HXX_ 