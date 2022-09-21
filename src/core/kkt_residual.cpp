#include "robotoc/core/kkt_residual.hpp"


namespace robotoc {

std::ostream& operator<<(std::ostream& os, const KKTResidual& kkt_residual) {
  os << "KKT residual:" << std::endl;
  kkt_residual.disp(os);
  return os;
}

} // namespace robotoc 