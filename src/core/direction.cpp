#include "robotoc/core/direction.hpp"


namespace robotoc {

std::ostream& operator<<(std::ostream& os, const Direction& d) {
  os << "Direction:" << std::endl;
  for (int i=0; i<d.size(); ++i) {
    os << "i: " << i << "\n";
    os << d[i] << "\n";
  }
  return os;
}

} // namespace robotoc 