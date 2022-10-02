#include "robotoc/core/kkt_matrix.hpp"


namespace robotoc {

std::ostream& operator<<(std::ostream& os, const KKTMatrix& kkt_matrix) {
  os << "KKTMatrix:" << "\n";
  for (int i=0; i<kkt_matrix.size(); ++i) {
    os << "i: " << i << "\n";
    os << kkt_matrix[i] << "\n";
  }
  return os;
}

} // namespace robotoc 