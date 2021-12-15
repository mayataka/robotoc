#include "robotoc/robot/impulse_status.hpp"


namespace robotoc {

void ImpulseStatus::disp(std::ostream& os) const {
  os << "impulse status:" << std::endl;
  os << "  impulse id: " << contact_status_.contactId() << std::endl;
  os << "  active impulses: [";
  for (int i=0; i<maxPointContacts()-1; ++i) {
    if (isImpulseActive(i)) {
      os << i << ", ";
    }
  }
  if (isImpulseActive(maxPointContacts()-1)) {
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
                         const ImpulseStatus& impulse_status) {
  impulse_status.disp(os);
  return os;
}

} // namespace robotoc 