#include "robotoc/hybrid/contact_sequence.hpp"


namespace robotoc {

void ContactSequence::disp(std::ostream& os) const {
  int impulse_index = 0;
  int lift_index = 0;
  os << "contact sequence:" << std::endl;
  for (int event_index=0; event_index<numDiscreteEvents(); ++event_index) {
    os << "  contact phase: " << event_index << std::endl;
    os << contactStatus(event_index) << std::endl;
    os << "  event index: " << event_index << ", type: ";
    if (eventType(event_index) == DiscreteEventType::Impulse) {
      os << "impulse, time: " << impulseTime(impulse_index) << std::endl;
      os << impulseStatus(impulse_index) << std::endl;
      ++impulse_index;
    }
    else {
      os << "lift, time: " << liftTime(lift_index) << std::endl;
      ++lift_index;
    }
  }
  os << "  contact phase: " << numDiscreteEvents() << std::endl;
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