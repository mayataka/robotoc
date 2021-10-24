#include "robotoc/impulse/impulse_split_direction.hpp"


namespace robotoc {

void ImpulseSplitDirection::disp(std::ostream& os) const {
  os << "impulse split diretion:" << std::endl;
  os << "  dq = " << dq().transpose() << std::endl;
  os << "  dv = " << dv().transpose() << std::endl;
  os << "  ddv = " << ddv().transpose() << std::endl;
  if (dimi_ > 0) {
    os << "  df = " << df().transpose() << std::endl;
  }
  os << "  dlmd = " << dlmd().transpose() << std::endl;
  os << "  dgmm = " << dgmm().transpose() << std::endl;
  os << "  dbeta = " << dbeta().transpose() << std::endl;
  if (dimi_ > 0) {
    os << "  dmu = " << dmu().transpose() << std::endl;
  }
}


std::ostream& operator<<(std::ostream& os, const ImpulseSplitDirection& d) {
  d.disp(os);
  return os;
}

} // namespace robotoc 