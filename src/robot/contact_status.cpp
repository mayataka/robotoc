#include "idocp/robot/contact_status.hpp"

#include <iostream>

namespace idocp {

void ContactStatus::showInfo() const {
  std::cout << "----- The contact status -----" << std::endl;
  std::cout << "active contacts: [";
  for (int i=0; i<maxPointContacts(); ++i) {
    if (isContactActive(i)) {
      std::cout << i << ", ";
    }
  }
  std::cout << "]" << std::endl;
  std::cout << "contact points: [";
  for (int i=0; i<maxPointContacts(); ++i) {
    std::cout << "[" << contactPoint(i).transpose() << "], ";
  }
  std::cout << "]" << std::endl;
  std::cout << "------------------------------" << std::endl;
}

} // namespace idocp 