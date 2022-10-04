#include "robotoc/core/solution.hpp"


namespace robotoc {

std::ostream& operator<<(std::ostream& os, const Solution& s) {
  os << "Solution:" << "\n";
  for (int i=0; i<s.size(); ++i) {
    os << "i: " << i << "\n";
    os << s[i] << "\n";
  }
  return os;
}

} // namespace robotoc 