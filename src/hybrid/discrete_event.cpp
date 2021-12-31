#include "robotoc/hybrid/discrete_event.hpp"


namespace robotoc {

void DiscreteEvent::disp(std::ostream& os) const {
  os << "event type: ";
  if (event_type_ == DiscreteEventType::Impulse) {
    os << "Impulse" << std::endl;
  }
  else if (event_type_ == DiscreteEventType::Lift) {
    os << "Lift" << std::endl;
  }
  else {
    os << "None" << std::endl;
  }
  os << "max num contacts: " << max_num_contacts_ << std::endl;
  os << "pre contact status:" << std::endl;
  os << pre_contact_status_ << std::endl;
  os << "post contact status:" << std::endl;
  os << post_contact_status_ << std::endl;
  os << "impulse status:" << std::endl;
  os << impulse_status_ << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const DiscreteEvent& discrete_event) {
  discrete_event.disp(os);
  return os;
}

} // namespace robotoc 