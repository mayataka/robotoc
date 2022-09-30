#include "robotoc/planner/discrete_event.hpp"


namespace robotoc {

void DiscreteEvent::disp(std::ostream& os) const {
  auto discreteEventTypeToString = [](const DiscreteEventType& type) {
    switch (type)
    {
    case DiscreteEventType::Impulse:
      return "Impulse";
      break;
    case DiscreteEventType::Lift:
      return "Impulse";
      break;
    default:
      return "None";
      break;
    }
  };
  os << "event type: " << discreteEventTypeToString(event_type_) << "\n";
  os << "max num contacts: " << max_num_contacts_ << "\n";
  os << "pre contact status:" << "\n";
  os << pre_contact_status_ << "\n";
  os << "post contact status:" << "\n";
  os << post_contact_status_ << "\n";
  os << "impulse status:" << "\n";
  os << impulse_status_ << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const DiscreteEvent& discrete_event) {
  discrete_event.disp(os);
  return os;
}

} // namespace robotoc 