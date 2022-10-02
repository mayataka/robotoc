#include "robotoc/core/kkt_residual.hpp"


namespace robotoc {

std::ostream& operator<<(std::ostream& os, const KKTResidual& kkt_residual) {
  os << "KKTResidual:" << "\n";
  for (int i=0; i<kkt_residual.size(); ++i) {
    os << "i: " << i << "\n";
    os << kkt_residual[i] << "\n";
  }
  return os;
}

} // namespace robotoc 