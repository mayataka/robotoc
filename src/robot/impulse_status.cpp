#include "robotoc/robot/impulse_status.hpp"


namespace robotoc {

void ImpulseStatus::disp(std::ostream& os) const {
  os << "impulse status:" << std::endl;
  os << "  impulse mode id: " << contact_status_.contactModeId() << std::endl;
  os << "  active impulses: [";
  for (int i=0; i<maxNumContacts()-1; ++i) {
    if (isImpulseActive(i)) {
      os << i << ", ";
    }
  }
  if (isImpulseActive(maxNumContacts()-1)) {
    os << maxNumContacts()-1;
  }
  os << "]" << std::endl;
  os << "  contact positions: [";
  for (int i=0; i<maxNumContacts()-1; ++i) {
    os << "[" << contactPosition(i).transpose() << "], ";
  }
  os << "[" << contactPosition(maxNumContacts()-1).transpose() << "]";
  os << "]" << std::endl;
  os << "  contact rotations: [";
  for (int i=0; i<maxNumContacts()-1; ++i) {
    os << "[" << contactRotation(i).row(0) << "]  ";
  }
  os << "[" << contactRotation(maxNumContacts()-1).row(0) << "]" << std::endl;
  os << "                               ";
  for (int i=0; i<maxNumContacts()-1; ++i) {
    os << "[" << contactRotation(i).row(1) << "]  ";
  }
  os << "[" << contactRotation(maxNumContacts()-1).row(1) << "]" << std::endl;
  os << "                               ";
  for (int i=0; i<maxNumContacts()-1; ++i) {
    os << "[" << contactRotation(i).row(2) << "], ";
  }
  os << "[" << contactRotation(maxNumContacts()-1).row(2) << "]";
  os << "]" << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const ImpulseStatus& impulse_status) {
  impulse_status.disp(os);
  return os;
}

} // namespace robotoc 