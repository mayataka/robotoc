#include "robotoc/robot/contact_status.hpp"


namespace robotoc {

void ContactStatus::disp(std::ostream& os) const {
  os << "contact status:" << std::endl;
  os << "  contact id: " << contact_id_ << std::endl;
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
  os << "]" << std::endl;
  os << "  contact surfaces rotations: [";
  for (int i=0; i<maxPointContacts()-1; ++i) {
    os << "[" << contactSurfaceRotation(i).row(0) << "]  ";
  }
  os << "[" << contactSurfaceRotation(maxPointContacts()-1).row(0) << "]" << std::endl;
  os << "                               ";
  for (int i=0; i<maxPointContacts()-1; ++i) {
    os << "[" << contactSurfaceRotation(i).row(1) << "]  ";
  }
  os << "[" << contactSurfaceRotation(maxPointContacts()-1).row(1) << "]" << std::endl;
  os << "                               ";
  for (int i=0; i<maxPointContacts()-1; ++i) {
    os << "[" << contactSurfaceRotation(i).row(2) << "], ";
  }
  os << "[" << contactSurfaceRotation(maxPointContacts()-1).row(2) << "]";
  os << "]" << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const ContactStatus& contact_status) {
  contact_status.disp(os);
  return os;
}

} // namespace robotoc 