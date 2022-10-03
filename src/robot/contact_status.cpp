#include "robotoc/robot/contact_status.hpp"


namespace robotoc {

void ContactStatus::disp(std::ostream& os) const {
  os << "ContactStatus:" << "\n";
  os << "  active contacts: [";
  for (int i=0; i<maxNumContacts()-1; ++i) {
    if (isContactActive(i)) {
      os << i << ", ";
    }
  }
  if (isContactActive(maxNumContacts()-1)) {
    os << maxNumContacts()-1;
  }
  os << "]" << "\n";
  os << "  contact positions: [";
  for (int i=0; i<maxNumContacts()-1; ++i) {
    os << "[" << contactPosition(i).transpose() << "], ";
  }
  os << "[" << contactPosition(maxNumContacts()-1).transpose() << "]";
  os << "]" << "\n";
  os << "  contact surfaces rotations: [";
  for (int i=0; i<maxNumContacts()-1; ++i) {
    os << "[" << contactRotation(i).row(0) << "]  ";
  }
  os << "[" << contactRotation(maxNumContacts()-1).row(0) << "]" << "\n";
  os << "                               ";
  for (int i=0; i<maxNumContacts()-1; ++i) {
    os << "[" << contactRotation(i).row(1) << "]  ";
  }
  os << "[" << contactRotation(maxNumContacts()-1).row(1) << "]" << "\n";
  os << "                               ";
  for (int i=0; i<maxNumContacts()-1; ++i) {
    os << "[" << contactRotation(i).row(2) << "], ";
  }
  os << "[" << contactRotation(maxNumContacts()-1).row(2) << "]";
  os << "]" << "\n";
  os << "  friction coefficients: [";
  for (int i=0; i<maxNumContacts()-1; ++i) {
    os << "[" << frictionCoefficient(i) << "], ";
  }
  os << "[" << frictionCoefficient(maxNumContacts()-1) << "]";
  os << "]" << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const ContactStatus& contact_status) {
  contact_status.disp(os);
  return os;
}

} // namespace robotoc 