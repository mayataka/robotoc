#include "robotoc/robot/contact_status.hpp"


namespace robotoc {

void ContactStatus::disp(std::ostream& os) const {
  os << "contact status:" << std::endl;
  os << "  active contacts: [";
  for (int i=0; i<maxPointContacts()-1; ++i) {
    if (isContactActive(i)) {
      os << i << ", ";
    }
  }
  if (isContactActive(maxPointContacts()-1)) {
    os << maxPointContacts()-1;
  }
  os << "]" << std::endl;
  os << "  contact points: [";
  for (int i=0; i<maxPointContacts()-1; ++i) {
    os << "[" << contactPoint(i).transpose() << "], ";
  }
  os << "[" << contactPoint(maxPointContacts()-1).transpose() << "]";
  os << "]" << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const ContactStatus& contact_status) {
  contact_status.disp(os);
  return os;
}

} // namespace robotoc 