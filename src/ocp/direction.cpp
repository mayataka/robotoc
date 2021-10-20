#include "idocp/ocp/direction.hpp"


namespace idocp {

std::ostream& operator<<(std::ostream& os, const Direction& d) {
  os << "direction:" << std::endl;
  d.disp(os);
  return os;
}

} // namespace idocp 