#ifndef IDOCP_CONTACT_SEQUENCE_HXX_
#define IDOCP_CONTACT_SEQUENCE_HXX_ 

#include "idocp/hybrid/contact_sequence.hpp"

#include <stdexcept>
#include <cassert>
#include <algorithm>

namespace idocp {

inline ContactSequence::ContactSequence(const Robot& robot, 
                                        const int max_num_events)
  : max_num_events_(max_num_events),
    num_impulse_events_(0), 
    num_lift_events_(0),
    contact_status_sequence_(max_num_events+1, robot.createContactStatus()),
    impulse_event_sequence_(max_num_events, DiscreteEvent(robot)),
    lift_event_sequence_(max_num_events, DiscreteEvent(robot)) {
  try {
    if (max_num_events <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline ContactSequence::ContactSequence()
  : max_num_events_(0),
    num_impulse_events_(0), 
    num_lift_events_(0),
    contact_status_sequence_(),
    impulse_event_sequence_(),
    lift_event_sequence_() {
}


inline ContactSequence::~ContactSequence() {
}


inline void ContactSequence::setContactStatusUniformly(
    const ContactStatus& contact_status) {
  for (auto& e : contact_status_sequence_) {
    e.set(contact_status);
  }
  for (auto& e : impulse_event_sequence_) {
    e.disableDiscreteEvent();
  }
  for (auto& e : lift_event_sequence_) {
    e.disableDiscreteEvent();
  }
  num_impulse_events_ = 0;
  num_lift_events_ = 0;
}


inline void ContactSequence::pushBackDiscreteEvent(
    const DiscreteEvent& discrete_event) {
  try {
    if (!discrete_event.existDiscreteEvent()) {
      throw std::runtime_error(
          "discrete_event.existDiscreteEvent() must be true!");
    }
    const int num_discrete_events = numImpulseEvents() + numLiftEvents();
    if (discrete_event.preContactStatus() 
          != contactStatus(num_discrete_events)) {
      throw std::runtime_error(
          "discrete_event.preContactStatus() is not consistent with this!");
    }
    if (num_discrete_events >= max_num_events_) {
      throw std::runtime_error(
          "Number of discrete events exceeds predefined max_num_events!");
    }
    if (numImpulseEvents() > 0 && numLiftEvents() > 0) {
      const double max_event_time = std::max(impulseTime(numImpulseEvents()-1),
                                             liftTime(numLiftEvents()-1));
      if (discrete_event.eventTime <= max_event_time) {
        throw std::runtime_error(
            "discrete_event.eventTime must be larger than max event time!");
      }
    }
    else if (numImpulseEvents() > 0 && numLiftEvents() == 0) {
      if (discrete_event.eventTime <= impulseTime(numImpulseEvents()-1)) {
        throw std::runtime_error(
            "discrete_event.eventTime must be larger than max event time!");
      }
    }
    else if (numImpulseEvents() == 0 && numLiftEvents() > 0) {
      if (discrete_event.eventTime <= liftTime(numLiftEvents()-1)) {
        throw std::runtime_error(
            "discrete_event.eventTime must be larger than max event time!");
      }
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  if (discrete_event.existImpulse()) {
    impulse_event_sequence_[num_impulse_events_] = discrete_event;
    ++num_impulse_events_;
  }
  else if (discrete_event.existLift()) {
    lift_event_sequence_[num_lift_events_] = discrete_event;
    ++num_lift_events_;
  }
  const int num_discrete_events = num_impulse_events_ + num_lift_events_;
  contact_status_sequence_[num_discrete_events].set(
      discrete_event.postContactStatus());
}


inline void ContactSequence::pushBackDiscreteEvent(
    const ContactStatus& contact_status, const double event_time) {
  const int num_discrete_events = num_impulse_events_ + num_lift_events_;
  DiscreteEvent discrete_event(contact_status_sequence_[num_discrete_events], 
                               contact_status);
  discrete_event.eventTime = event_time;
  pushBackDiscreteEvent(discrete_event);
}


inline void ContactSequence::popBackDiscreteEvent() {
  if (num_impulse_events_ > 0 && num_lift_events_ > 0) {
    const double max_impulse_time 
        = impulse_event_sequence_[num_impulse_events_-1].eventTime;
    const double max_lift_time 
        = lift_event_sequence_[num_lift_events_-1].eventTime;
    if (max_impulse_time > max_lift_time) {
      popBackImpulseEvent();
    }
    else {
      popBackLiftEvent();
    }
  }
  else if (num_impulse_events_ > 0 && num_lift_events_ == 0) {
    popBackImpulseEvent();
  }
  else if (num_impulse_events_ == 0 && num_lift_events_ > 0) {
    popBackLiftEvent();
  }
}


inline void ContactSequence::popFrontDiscreteEvent() {
  const int num_discrete_events = num_impulse_events_ + num_lift_events_;
  if (num_discrete_events > 0) {
    for (int i=0; i<=num_discrete_events-1; ++i) {
      contact_status_sequence_[i] = contact_status_sequence_[i+1];
    }
    if (num_impulse_events_ > 0 && num_lift_events_ > 0) {
      const double min_impulse_time = impulseTime(0);
      const double min_lift_time = liftTime(0);
      if (min_impulse_time < min_lift_time) {
        popFrontImpulseEvent();
      }
      else {
        popFrontLiftEvent();
      }
    }
    else if (num_impulse_events_ > 0 && num_lift_events_ == 0) {
      popFrontImpulseEvent();
    }
    else if (num_impulse_events_ == 0 && num_lift_events_ > 0) {
      popFrontLiftEvent();
    }
  }
}


inline void ContactSequence::shiftImpulseEvent(const int impulse_index, 
                                               const double event_time) {
  try {
    if (num_impulse_events_ <= 0) {
      throw std::runtime_error(
          "numImpulseEvents() must be positive when calling this method!");
    }
    if (impulse_index < 0) {
      throw std::runtime_error("impulse_index must be non-negative!");
    }
    if (impulse_index >= num_impulse_events_) {
      throw std::runtime_error(
          "impulse_index must be less than numImpulseEvents()!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  impulse_event_sequence_[impulse_index] = event_time;
}


inline void ContactSequence::shiftLiftEvent(const int lift_index, 
                                            const double event_time) {
  try {
    if (num_lift_events_ <= 0) {
      throw std::runtime_error(
          "numLiftEvents() must be positive when calling this method!");
    }
    if (lift_index < 0) {
      throw std::runtime_error("lift_index must be non-negative!");
    }
    if (lift_index >= num_lift_events_) {
      throw std::runtime_error("lift_index must be less than numLiftEvents()!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  lift_event_sequence_[lift_index] = event_time;
}


inline void ContactSequence::setContactPointsToImpulseEvent(
    const int impulse_index, Robot& robot, const Eigen::VectorXd& q) {
  try {
    if (impulse_index >= num_impulse_events_) {
      throw std::runtime_error(
          "impulse_index must be less than numImpulseEvents()!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  impulse_event_sequence_[impulse_index].setContactPoints(robot);
}


inline void ContactSequence::setContactPointsToLiftEvent(
    const int lift_index, Robot& robot, const Eigen::VectorXd& q) {
  try {
    if (lift_index >= num_lift_events_) {
      throw std::runtime_error(
          "lift_index must be less than numLiftEvents()!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  lift_event_sequence_[lift_index].setContactPoints(robot);
}


inline void ContactSequence::setContactPointsToContatStatus(
    const int contact_phase, Robot& robot, const Eigen::VectorXd& q) {
  try {
    if (contact_phase > num_impulse_events_+num_lift_events_) {
      throw std::runtime_error(
          "contact_phase must be less than numImpulseEvents()+numLiftEvents()!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  robot.setContactPoints(contact_status_sequence_[contact_phase]);
}


inline void ContactSequence::setContactPointsUniformly(
    Robot& robot, const Eigen::VectorXd& q) {
  for (int i=0; i<=num_impulse_events_+num_lift_events_; ++i) {
    robot.setContactPoints(contact_status_sequence_[i]);
  }
  for (int i=0; i<num_impulse_events_; ++i) {
    impulse_event_sequence_[i].setContactPoints(robot);
  }
  for (int i=0; i<num_lift_events_; ++i) {
    lift_event_sequence_[i].setContactPoints(robot);
  }
}


inline int ContactSequence::numImpulseEvents() const {
  return num_impulse_events_;
}


inline int ContactSequence::numLiftEvents() const {
  return num_lift_events_;
}


inline const ContactStatus& ContactSequence::contactStatus(
    const int contact_phase) const {
  assert(contact_phase >= 0);
  assert(contact_phase <= num_impulse_events_+num_lift_events_);
  return contact_status_sequence_[contact_phase];
}


inline const ImpulseStatus& ContactSequence::impulseStatus(
    const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < num_impulse_events_);
  assert(impulse_event_sequence_[impulse_index].existImpulse());
  return impulse_event_sequence_[impulse_index].impulseStatus();
}


inline double ContactSequence::impulseTime(const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < num_impulse_events_);
  return impulse_event_sequence_[impulse_index].eventTime;
}


inline double ContactSequence::liftTime(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < num_lift_events_);
  return lift_event_sequence_[lift_index].eventTime;
}


inline void ContactSequence::popBackImpulseEvent() {
  if (num_impulse_events_ > 0)  {
    --num_impulse_events_;
  }
}


inline void ContactSequence::popBackLiftEvent() {
  if (num_lift_events_ > 0)  {
    --num_lift_events_;
  }
}


inline void ContactSequence::popFrontImpulseEvent() {
  if (num_impulse_events_ > 0)  {
    for (int i=0; i<=num_impulse_events_-2; ++i) {
      impulse_event_sequence_[i] = impulse_event_sequence_[i+1];
    }
    --num_impulse_events_;
  }
}


inline void ContactSequence::popFrontLiftEvent() {
  if (num_lift_events_ > 0)  {
    for (int i=0; i<=num_lift_events_-2; ++i) {
      lift_event_sequence_[i] = lift_event_sequence_[i+1];
    }
    --num_lift_events_;
  }
}

} // namespace idocp 

#endif // IDOCP_CONTACT_SEQUENCE_HXX_ 