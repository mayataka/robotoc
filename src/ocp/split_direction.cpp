#include "idocp/ocp/split_direction.hpp"


namespace idocp {

void SplitDirection::disp(std::ostream& os) const {
  os << "split diretion:" << std::endl;
  os << "  dq = " << dq().transpose() << std::endl;
  os << "  dv = " << dv().transpose() << std::endl;
  os << "  du = " << du.transpose() << std::endl;
  os << "  da = " << da().transpose() << std::endl;
  if (dimf_ > 0) {
    os << "  df = " << df().transpose() << std::endl;
  }
  os << "  dlmd = " << dlmd().transpose() << std::endl;
  os << "  dgmm = " << dgmm().transpose() << std::endl;
  os << "  dbeta = " << dbeta().transpose() << std::endl;
  if (dimf_ > 0) {
    os << "  dmu = " << dmu().transpose() << std::endl;
  }
  if (dimi_ > 0) {
    os << "  dxi = " << dxi().transpose() << std::flush;
  }
}


std::ostream& operator<<(std::ostream& os, const SplitDirection& d) {
  d.disp(os);
  return os;
}

} // namespace idocp 