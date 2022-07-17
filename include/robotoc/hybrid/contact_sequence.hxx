#ifndef ROBOTOC_CONTACT_SEQUENCE_HXX_
#define ROBOTOC_CONTACT_SEQUENCE_HXX_ 

#include "robotoc/hybrid/contact_sequence.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>
#include <algorithm>

namespace robotoc {

inline ContactSequence::ContactSequence(const Robot& robot, 
                                        const int reserved_num_discrete_events)
  : reserved_num_discrete_events_(reserved_num_discrete_events),
    default_contact_status_(robot.createContactStatus()),
    contact_statuses_(2*reserved_num_discrete_events+1),
    impulse_events_(reserved_num_discrete_events),
    event_index_impulse_(reserved_num_discrete_events), 
    event_index_lift_(reserved_num_discrete_events),
    event_time_(2*reserved_num_discrete_events),
    impulse_time_(reserved_num_discrete_events),
    lift_time_(reserved_num_discrete_events),
    is_impulse_event_(2*reserved_num_discrete_events),
    sto_impulse_(reserved_num_discrete_events), 
    sto_lift_(reserved_num_discrete_events) {
  try {
    if (reserved_num_discrete_events < 0) {
      throw std::out_of_range("invalid argument: reserved_num_discrete_events must be non-negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  clear_all();
  contact_statuses_.push_back(default_contact_status_);
}


inline ContactSequence::ContactSequence()
  : reserved_num_discrete_events_(0),
    default_contact_status_(),
    contact_statuses_(),
    impulse_events_(),
    event_index_impulse_(), 
    event_index_lift_(),
    event_time_(),
    impulse_time_(),
    lift_time_(),
    is_impulse_event_(),
    sto_impulse_(),
    sto_lift_() {
}


inline ContactSequence::~ContactSequence() {
}


inline std::shared_ptr<ContactSequence> ContactSequence::clone() const {
  return std::make_shared<ContactSequence>(*this); 
} 


inline void ContactSequence::init(
    const ContactStatus& contact_status) {
  clear_all();
  contact_statuses_.push_back(contact_status);
}


inline void ContactSequence::push_back(const DiscreteEvent& discrete_event, 
                                       const double event_time,
                                       const bool sto) {
  try {
    if (numContactPhases() == 0) {
      throw std::runtime_error(
          "Call init() before calling push_back()!");
    }
    if (!discrete_event.existDiscreteEvent()) {
      throw std::runtime_error(
          "discrete_event.existDiscreteEvent() must be true!");
    }
    if (discrete_event.preContactStatus() != contact_statuses_.back()) {
      throw std::runtime_error(
          "discrete_event.preContactStatus() is not consistent with the last contact status!");
    }
    if (numImpulseEvents() > 0 || numLiftEvents() > 0) {
      if (event_time <= event_time_.back()) {
        throw std::runtime_error(
            "The input event_time " + std::to_string(event_time) 
            + " must be larger than the last event time=" 
            + std::to_string(event_time_.back()) + "!");
      }
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::cerr << "c.f. the current contact sequence is " << "\n";
    std::cerr << *this << "\n";
    std::exit(EXIT_FAILURE);
  }
  contact_statuses_.push_back(discrete_event.postContactStatus());
  event_time_.push_back(event_time);
  if (discrete_event.existImpulse()) {
    impulse_events_.push_back(discrete_event);
    event_index_impulse_.push_back(contact_statuses_.size()-2);
    impulse_time_.push_back(event_time);
    is_impulse_event_.push_back(true);
    sto_impulse_.push_back(sto);
  }
  else {
    event_index_lift_.push_back(contact_statuses_.size()-2);
    lift_time_.push_back(event_time);
    is_impulse_event_.push_back(false);
    sto_lift_.push_back(sto);
  }
  if (reserved_num_discrete_events_ < numDiscreteEvents()) {
    reserved_num_discrete_events_ = numDiscreteEvents();
  }
}


inline void ContactSequence::push_back(const ContactStatus& contact_status, 
                                       const double switching_time,
                                       const bool sto) {
  DiscreteEvent discrete_event(contact_statuses_.back(), contact_status);
  push_back(discrete_event, switching_time, sto);
}


inline void ContactSequence::pop_back() {
  if (numDiscreteEvents() > 0) {
    if (is_impulse_event_.back()) {
      impulse_events_.pop_back();
      event_index_impulse_.pop_back();
      impulse_time_.pop_back();
      sto_impulse_.pop_back();
    }
    else {
      event_index_lift_.pop_back();
      lift_time_.pop_back();
      sto_lift_.pop_back();
    }
    event_time_.pop_back();
    is_impulse_event_.pop_back();
    contact_statuses_.pop_back();
  }
  else if (numContactPhases() > 0) {
    assert(numContactPhases() == 1);
    contact_statuses_.pop_back();
    contact_statuses_.push_back(default_contact_status_);
  }
}


inline void ContactSequence::pop_front() {
  if (numDiscreteEvents() > 0) {
    if (is_impulse_event_.front()) {
      impulse_events_.pop_front();
      event_index_impulse_.pop_front();
      impulse_time_.pop_front();
      sto_impulse_.pop_front();
    }
    else {
      event_index_lift_.pop_front();
      lift_time_.pop_front();
      sto_lift_.pop_front();
    }
    event_time_.pop_front();
    is_impulse_event_.pop_front();
    contact_statuses_.pop_front();
    for (auto& e : event_index_impulse_) { e -= 1; }
    for (auto& e : event_index_lift_) { e -= 1; }
  }
  else if (numContactPhases() > 0) {
    assert(numContactPhases() == 1);
    contact_statuses_.pop_front();
    contact_statuses_.push_back(default_contact_status_);
  }
}


inline void ContactSequence::setImpulseTime(const int impulse_index, 
                                            const double impulse_time) {
  try {
    if (numImpulseEvents() <= 0) {
      throw std::runtime_error(
          "numImpulseEvents() must be positive when calling this method!");
    }
    if (impulse_index < 0) {
      throw std::runtime_error("impulse_index must be non-negative!");
    }
    if (impulse_index >= numImpulseEvents()) {
      throw std::runtime_error(
          "The input impulse_index " + std::to_string(impulse_index) 
          + " must be less than numImpulseEvents()=" 
          + std::to_string(numImpulseEvents()) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::cerr << "c.f. the current contact sequence is " << "\n";
    std::cerr << *this << "\n";
    std::exit(EXIT_FAILURE);
  }
  impulse_time_[impulse_index] = impulse_time;
  event_time_[event_index_impulse_[impulse_index]] = impulse_time;
}


inline void ContactSequence::setLiftTime(const int lift_index, 
                                         const double lift_time) {
  try {
    if (numLiftEvents() <= 0) {
      throw std::runtime_error(
          "numLiftEvents() must be positive when calling this method!");
    }
    if (lift_index < 0) {
      throw std::runtime_error("lift_index must be non-negative!");
    }
    if (lift_index >= numLiftEvents()) {
      throw std::runtime_error(
          "The input lift_index " + std::to_string(lift_index) 
          + " must be less than numLiftEvents()=" 
          + std::to_string(numLiftEvents()) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::cerr << "c.f. the current contact sequence is " << "\n";
    std::cerr << *this << "\n";
    std::exit(EXIT_FAILURE);
  }
  lift_time_[lift_index] = lift_time;
  event_time_[event_index_lift_[lift_index]] = lift_time;
}


inline bool ContactSequence::isEventTimeConsistent() const {
  bool is_consistent = true;
  if (numDiscreteEvents() > 0) {
    for (int event_index=1; event_index<numDiscreteEvents(); ++event_index) {
      if (event_time_[event_index] <= event_time_[event_index-1]) {
            "event_time[" + std::to_string(event_index) + "]=" 
            + std::to_string(event_time_[event_index]) 
            + " must be larger than event_time_[" 
            + std::to_string(event_index-1) + "]=" 
            + std::to_string(event_time_[event_index-1]) + "!";
        is_consistent = false;
      }
    }
  }
  return is_consistent;
}


inline void ContactSequence::setContactPlacements(
    const int contact_phase, 
    const std::vector<Eigen::Vector3d>& contact_positions) {
  try {
    if (contact_phase >= numContactPhases()) {
      throw std::runtime_error(
          "The input contact_phase " + std::to_string(contact_phase) 
          + " must be smaller than numContactPhases()" 
          + std::to_string(numContactPhases()) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::cerr << "c.f. the current contact sequence is " << "\n";
    std::cerr << *this << "\n";
    std::exit(EXIT_FAILURE);
  }
  contact_statuses_[contact_phase].setContactPlacements(contact_positions);
  if (contact_phase > 0) {
    if (is_impulse_event_[contact_phase-1]) {
      for (int impulse_index=0; ; ++impulse_index) {
        assert(impulse_index < numImpulseEvents());
        if (event_index_impulse_[impulse_index] == contact_phase-1) {
          impulse_events_[impulse_index].setContactPlacements(contact_positions);
          break;
        }
      }
    }
  }
}


inline void ContactSequence::setContactPlacements(
    const int contact_phase, 
    const std::vector<Eigen::Vector3d>& contact_positions,
    const std::vector<Eigen::Matrix3d>& contact_rotations) {
  try {
    if (contact_phase >= numContactPhases()) {
      throw std::runtime_error(
          "The input contact_phase " + std::to_string(contact_phase) 
          + " must be smaller than numContactPhases()" 
          + std::to_string(numContactPhases()) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::cerr << "c.f. the current contact sequence is " << "\n";
    std::cerr << *this << "\n";
    std::exit(EXIT_FAILURE);
  }
  contact_statuses_[contact_phase].setContactPlacements(contact_positions, 
                                                        contact_rotations);
  if (contact_phase > 0) {
    if (is_impulse_event_[contact_phase-1]) {
      for (int impulse_index=0; ; ++impulse_index) {
        assert(impulse_index < numImpulseEvents());
        if (event_index_impulse_[impulse_index] == contact_phase-1) {
          impulse_events_[impulse_index].setContactPlacements(contact_positions,
                                                              contact_rotations);
          break;
        }
      }
    }
  }
}


inline void ContactSequence::setContactPlacements(
    const int contact_phase, const aligned_vector<SE3>& contact_placements) {
  try {
    if (contact_phase >= numContactPhases()) {
      throw std::runtime_error(
          "The input contact_phase " + std::to_string(contact_phase) 
          + " must be smaller than numContactPhases()" 
          + std::to_string(numContactPhases()) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::cerr << "c.f. the current contact sequence is " << "\n";
    std::cerr << *this << "\n";
    std::exit(EXIT_FAILURE);
  }
  contact_statuses_[contact_phase].setContactPlacements(contact_placements);
  if (contact_phase > 0) {
    if (is_impulse_event_[contact_phase-1]) {
      for (int impulse_index=0; ; ++impulse_index) {
        assert(impulse_index < numImpulseEvents());
        if (event_index_impulse_[impulse_index] == contact_phase-1) {
          impulse_events_[impulse_index].setContactPlacements(contact_placements);
          break;
        }
      }
    }
  }
}


inline int ContactSequence::numImpulseEvents() const {
  return impulse_events_.size();
}


inline int ContactSequence::numLiftEvents() const {
  return (numDiscreteEvents()-numImpulseEvents());
}


inline int ContactSequence::numDiscreteEvents() const {
  return (numContactPhases()-1);
}


inline int ContactSequence::numContactPhases() const {
  return contact_statuses_.size();
}


inline const ContactStatus& ContactSequence::contactStatus(
    const int contact_phase) const {
  assert(contact_phase >= 0);
  assert(contact_phase < numContactPhases());
  return contact_statuses_[contact_phase];
}


inline const ImpulseStatus& ContactSequence::impulseStatus(
    const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < numImpulseEvents());
  return impulse_events_[impulse_index].impulseStatus();
}


inline double ContactSequence::impulseTime(const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < numImpulseEvents());
  return impulse_time_[impulse_index];
}


inline double ContactSequence::liftTime(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < numLiftEvents());
  return lift_time_[lift_index];
}


inline bool ContactSequence::isSTOEnabledImpulse(const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < numImpulseEvents());
  return sto_impulse_[impulse_index];
}


inline bool ContactSequence::isSTOEnabledLift(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < numLiftEvents());
  return sto_lift_[lift_index];
}


inline DiscreteEventType ContactSequence::eventType(
    const int event_index) const {
  assert(event_index >= 0);
  assert(event_index < numImpulseEvents()+numLiftEvents());
  if (is_impulse_event_[event_index]) return DiscreteEventType::Impulse;
  else return DiscreteEventType::Lift;
}


inline const std::deque<double>& ContactSequence::eventTimes() const {
  return event_time_;
}


inline void ContactSequence::reserve(const int reserved_num_discrete_events) {
  if (reserved_num_discrete_events_ < reserved_num_discrete_events) {
    reserveDeque(contact_statuses_, 2*reserved_num_discrete_events+1);
    reserveDeque(impulse_events_, reserved_num_discrete_events);
    reserveDeque(event_index_impulse_, reserved_num_discrete_events);
    reserveDeque(event_index_lift_, reserved_num_discrete_events);
    reserveDeque(event_time_, 2*reserved_num_discrete_events);
    reserveDeque(impulse_time_, reserved_num_discrete_events);
    reserveDeque(lift_time_, reserved_num_discrete_events);
    reserveDeque(is_impulse_event_, 2*reserved_num_discrete_events);
    reserveDeque(sto_impulse_, reserved_num_discrete_events);
    reserveDeque(sto_lift_, reserved_num_discrete_events);
    reserved_num_discrete_events_ = reserved_num_discrete_events;
  }
}


inline int ContactSequence::reservedNumDiscreteEvents() const {
  return reserved_num_discrete_events_;
}


inline void ContactSequence::clear_all() {
  contact_statuses_.clear();
  impulse_events_.clear();
  event_index_impulse_.clear(); 
  event_index_lift_.clear();
  event_time_.clear();
  impulse_time_.clear();
  lift_time_.clear();
  is_impulse_event_.clear();
  sto_impulse_.clear();
  sto_lift_.clear();
}

} // namespace robotoc 

#endif // ROBOTOC_CONTACT_SEQUENCE_HXX_ 