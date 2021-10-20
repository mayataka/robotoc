#include "idocp/ocp/solution.hpp"


namespace idocp {

std::ostream& operator<<(std::ostream& os, const Solution& s) {
  os << "solution:" << std::endl;
  s.disp(os);
  return os;
}

} // namespace idocp 