#include "idocp/robot/impulse_status.hpp"

#include <iostream>

namespace idocp {

void ImpulseStatus::showInfo() const {
  std::cout << "----- The impulse status -----" << std::endl;
  std::cout << "active impulses: [";
  for (int i=0; i<maxPointContacts(); ++i) {
    if (isImpulseActive(i)) {
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