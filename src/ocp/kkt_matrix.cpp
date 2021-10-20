#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

std::ostream& operator<<(std::ostream& os, const KKTMatrix& kkt_matrix) {
  os << "KKT matrix:" << std::endl;
  kkt_matrix.disp(os);
  return os;
}

} // namespace idocp 