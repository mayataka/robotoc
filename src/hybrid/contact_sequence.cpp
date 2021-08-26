#include "idocp/hybrid/contact_sequence.hpp"

#include <iostream>

namespace idocp {

void ContactSequence::showInfo() const {
  int impulse_index = 0;
  int lift_index = 0;
  std::cout << "----- The contact sequence -----" << std::endl;
  for (int event_index=0; event_index<numDiscreteEvents(); ++event_index) {
    std::cout << "contact phase: " << event_index << std::endl;
    contactStatus(event_index).showInfo();
    std::cout << "event index: " << event_index << ", type: ";
    if (eventType(event_index) == DiscreteEventType::Impulse) {
      std::cout << "impulse, time: " << impulseTime(impulse_index) << std::endl;
      impulseStatus(impulse_index).showInfo();
      ++impulse_index;
    }
    else {
      std::cout << "lift, time: " << liftTime(lift_index) << std::endl;
      ++lift_index;
    }
  }
  std::cout << "contact phase: " << numDiscreteEvents() << std::endl;
  contactStatus(numDiscreteEvents()).showInfo();
  std::cout << "--------------------------------" << std::endl;
}

} // namespace idocp 