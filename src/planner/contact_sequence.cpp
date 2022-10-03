#include "robotoc/planner/contact_sequence.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>
#include <algorithm>


namespace robotoc {

ContactSequence::ContactSequence(const Robot& robot, 
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
  if (reserved_num_discrete_events < 0) {
    throw std::out_of_range("[ContactSequence] invalid argument: reserved_num_discrete_events must be non-negative!");
  }
  clear();
  contact_statuses_.push_back(default_contact_status_);
}


ContactSequence::ContactSequence()
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


void ContactSequence::init(const ContactStatus& contact_status) {
  clear();
  contact_statuses_.push_back(contact_status);
}


void ContactSequence::push_back(const DiscreteEvent& discrete_event, 
                                const double event_time, const bool sto) {
  if (numContactPhases() == 0) {
    throw std::runtime_error(
        "[ContactSequence] call init() before calling push_back()!");
  }
  if (!discrete_event.existDiscreteEvent()) {
    throw std::runtime_error(
        "[ContactSequence] discrete_event.existDiscreteEvent() must be true!");
  }
  if (discrete_event.preContactStatus() != contact_statuses_.back()) {
    throw std::runtime_error(
        "[ContactSequence] discrete_event.preContactStatus() is not consistent with the last contact status!");
  }
  if (numImpulseEvents() > 0 || numLiftEvents() > 0) {
    if (event_time <= event_time_.back()) {
      throw std::runtime_error(
          "[ContactSequence] input event_time (" + std::to_string(event_time) 
          + ") must be larger than the last event time (" 
          + std::to_string(event_time_.back()) + ") !");
    }
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


void ContactSequence::push_back(const ContactStatus& contact_status, 
                                const double switching_time, const bool sto) {
  DiscreteEvent discrete_event(contact_statuses_.back(), contact_status);
  push_back(discrete_event, switching_time, sto);
}


void ContactSequence::pop_back() {
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


void ContactSequence::pop_front() {
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


void ContactSequence::setImpulseTime(const int impulse_index, 
                                     const double impulse_time) {
  if (numImpulseEvents() <= 0) {
    throw std::runtime_error(
        "[ContactSequence] numImpulseEvents() must be positive when calling setImpulseTime()!");
  }
  if (impulse_index < 0) {
    throw std::runtime_error("[ContactSequence] 'impulse_index' must be non-negative!");
  }
  if (impulse_index >= numImpulseEvents()) {
    throw std::runtime_error(
        "[ContactSequence] input 'impulse_index' (" + std::to_string(impulse_index) 
        + ") must be less than numImpulseEvents() (" 
        + std::to_string(numImpulseEvents()) + ") !");
  }
  impulse_time_[impulse_index] = impulse_time;
  event_time_[event_index_impulse_[impulse_index]] = impulse_time;
}


void ContactSequence::setLiftTime(const int lift_index, const double lift_time) {
  if (numLiftEvents() <= 0) {
    throw std::runtime_error(
        "[ContactSequence] numLiftEvents() must be positive when calling setLiftTime()!");
  }
  if (lift_index < 0) {
    throw std::runtime_error("[ContactSequence] 'lift_index' must be non-negative!");
  }
  if (lift_index >= numLiftEvents()) {
    throw std::runtime_error(
        "[ContactSequence] input 'lift_index' (" + std::to_string(lift_index) 
        + ") must be less than numLiftEvents() (" 
        + std::to_string(numLiftEvents()) + ") !");
  }
  lift_time_[lift_index] = lift_time;
  event_time_[event_index_lift_[lift_index]] = lift_time;
}


bool ContactSequence::isSTOEnabledImpulse(const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < numImpulseEvents());
  return sto_impulse_[impulse_index];
}


bool ContactSequence::isSTOEnabledLift(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < numLiftEvents());
  return sto_lift_[lift_index];
}


bool ContactSequence::isEventTimeConsistent() const {
  bool is_consistent = true;
  if (numDiscreteEvents() > 0) {
    for (int event_index=1; event_index<numDiscreteEvents(); ++event_index) {
      if (event_time_[event_index] <= event_time_[event_index-1]) {
        std::cerr << "[ContactSequence] event_time[" + std::to_string(event_index) + "] (" 
                        + std::to_string(event_time_[event_index]) 
                        + ") must be larger than event_time_[" 
                        + std::to_string(event_index-1) + "] (" 
                        + std::to_string(event_time_[event_index-1]) + ") !" << std::endl;
        is_consistent = false;
      }
    }
  }
  return is_consistent;
}


void ContactSequence::setContactPlacements(
    const int contact_phase, 
    const std::vector<Eigen::Vector3d>& contact_positions) {
  if (contact_phase >= numContactPhases()) {
    throw std::runtime_error(
        "[ContactSequence] input 'contact_phase' (" + std::to_string(contact_phase) 
        + ") must be smaller than numContactPhases() (" 
        + std::to_string(numContactPhases()) + ") !");
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


void ContactSequence::setContactPlacements(
    const int contact_phase, 
    const std::vector<Eigen::Vector3d>& contact_positions,
    const std::vector<Eigen::Matrix3d>& contact_rotations) {
  if (contact_phase >= numContactPhases()) {
    throw std::runtime_error(
        "[ContactSequence] input 'contact_phase' (" + std::to_string(contact_phase) 
        + ") must be smaller than numContactPhases() (" 
        + std::to_string(numContactPhases()) + ") !");
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


void ContactSequence::setContactPlacements(
    const int contact_phase, const aligned_vector<SE3>& contact_placements) {
  if (contact_phase >= numContactPhases()) {
    throw std::runtime_error(
        "[ContactSequence] input 'contact_phase' (" + std::to_string(contact_phase) 
        + ") must be smaller than numContactPhases() (" 
        + std::to_string(numContactPhases()) + ") !");
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


void ContactSequence::setFrictionCoefficients(
    const int contact_phase, const std::vector<double>& friction_coefficient) {
  if (contact_phase >= numContactPhases()) {
    throw std::runtime_error(
        "[ContactSequence] input 'contact_phase' (" + std::to_string(contact_phase) 
        + ") must be smaller than numContactPhases() (" 
        + std::to_string(numContactPhases()) + ") !");
  }
  contact_statuses_[contact_phase].setFrictionCoefficients(friction_coefficient);
  if (contact_phase > 0) {
    if (is_impulse_event_[contact_phase-1]) {
      for (int impulse_index=0; ; ++impulse_index) {
        assert(impulse_index < numImpulseEvents());
        if (event_index_impulse_[impulse_index] == contact_phase-1) {
          impulse_events_[impulse_index].setFrictionCoefficients(friction_coefficient);
          break;
        }
      }
    }
  }
}


template <typename T> 
void reserveDeque(std::deque<T>& deq, const int size) {
  assert(size >= 0);
  if (deq.empty()) {
    deq.resize(size);
  }
  else {
    const int current_size = deq.size();
    if (current_size < size) {
      while (deq.size() < size) {
        deq.push_back(deq.back());
      }
      while (deq.size() > current_size) {
        deq.pop_back();
      }
    }
  }
}


void ContactSequence::reserve(const int reserved_num_discrete_events) {
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


int ContactSequence::reservedNumDiscreteEvents() const {
  return reserved_num_discrete_events_;
}


void ContactSequence::clear() {
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



void ContactSequence::disp(std::ostream& os) const {
  int impulse_index = 0;
  int lift_index = 0;
  os << "contact sequence:" << "\n";
  for (int event_index=0; event_index<numDiscreteEvents(); ++event_index) {
    os << "  contact phase: " << event_index << "\n";
    os << contactStatus(event_index) << "\n";
    os << "  event index: " << event_index << ", type: ";
    if (eventType(event_index) == DiscreteEventType::Impact) {
      os << "impulse, time: " << impulseTime(impulse_index) 
         << ", sto: " << std::boolalpha << isSTOEnabledImpulse(impulse_index) <<  "\n";
      os << impulseStatus(impulse_index) << "\n";
      ++impulse_index;
    }
    else {
      os << "lift, time: " << liftTime(lift_index) 
         << ", sto: " << std::boolalpha << isSTOEnabledLift(lift_index) << "\n";
      ++lift_index;
    }
  }
  os << "  contact phase: " << numDiscreteEvents() << "\n";
  os << contactStatus(numDiscreteEvents()) << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const ContactSequence& contact_sequence) {
  contact_sequence.disp(os);
  return os;
}


std::ostream& operator<<(
    std::ostream& os, 
    const std::shared_ptr<ContactSequence>& contact_sequence) {
  contact_sequence->disp(os);
  return os;
}

} // namespace robotoc 