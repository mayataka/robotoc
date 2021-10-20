#include "idocp/impulse/impulse_split_kkt_matrix.hpp"


namespace idocp {

void ImpulseSplitKKTMatrix::disp(std::ostream& os) const {
  os << "impulse split KKT matrix:" << std::endl;
  os << "  Fxx = " << Fxx << std::endl;
  os << "  Fqq_prev = " << Fqq_prev << std::endl;
  os << "  Qxx = " << Qxx << std::endl;
  os << "  Qdvdv = " << Qdvdv << std::endl;
  if (dimi_ > 0) {
    os << "  Qff = " << Qff() << std::endl;
    os << "  Qqf = " << Qqf() << std::endl;
  }
}


std::ostream& operator<<(std::ostream& os, 
                         const ImpulseSplitKKTMatrix& kkt_matrix) {
  kkt_matrix.disp(os);
  return os;
}

} // namespace idocp 