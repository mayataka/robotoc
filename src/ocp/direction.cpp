#include "robotoc/ocp/direction.hpp"


namespace robotoc {

std::ostream& operator<<(std::ostream& os, const Direction& d) {
  os << "direction:" << std::endl;
  d.disp(os);
  return os;
}

} // namespace robotoc 