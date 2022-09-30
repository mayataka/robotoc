#include "robotoc/planner/contact_sequence.hpp"


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
  clear_all();
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


ContactSequence::~ContactSequence() {
}


void ContactSequence::disp(std::ostream& os) const {
  int impulse_index = 0;
  int lift_index = 0;
  os << "contact sequence:" << "\n";
  for (int event_index=0; event_index<numDiscreteEvents(); ++event_index) {
    os << "  contact phase: " << event_index << "\n";
    os << contactStatus(event_index) << "\n";
    os << "  event index: " << event_index << ", type: ";
    if (eventType(event_index) == DiscreteEventType::Impulse) {
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