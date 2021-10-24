#include "robotoc/ocp/solution.hpp"


namespace robotoc {

std::ostream& operator<<(std::ostream& os, const Solution& s) {
  os << "solution:" << std::endl;
  s.disp(os);
  return os;
}

} // namespace robotoc 