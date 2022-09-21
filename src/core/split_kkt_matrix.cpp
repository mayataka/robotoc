#include "robotoc/core/split_kkt_matrix.hpp"


namespace robotoc {

void SplitKKTMatrix::disp(std::ostream& os) const {
  os << "split KKT matrix:" << std::endl;
  os << "  Fxx = " << Fxx << std::endl;
  os << "  Fvu = " << Fvu << std::endl;
  os << "  Fqq_prev = " << Fqq_prev << std::endl;
  os << "  Qxx = " << Qxx << std::endl;
  os << "  Qxu = " << Qxu << std::endl;
  os << "  Quu = " << Quu << std::endl;
  os << "  Qaa = " << Qaa << std::endl;
  if (dimf_ > 0) {
    os << "  Qff = " << Qff() << std::endl;
    os << "  Qqf = " << Qqf() << std::endl;
  }
  os << "  fq = " << fq().transpose() << std::endl;
  os << "  fv = " << fv().transpose() << std::endl;
  os << "  Qtt = " << Qtt << std::endl;
  os << "  Qtt_prev = " << Qtt_prev << std::endl;
  os << "  hq = " << hq().transpose() << std::endl;
  os << "  hv = " << hv().transpose() << std::endl;
  os << "  hu = " << hu.transpose() << std::endl;
  os << "  ha = " << ha.transpose() << std::flush;
  if (dimf_ > 0) {
    os << std::endl;
    os << "  hf = " << hf().transpose() << std::flush;
  }
}


std::ostream& operator<<(std::ostream& os, const SplitKKTMatrix& kkt_matrix) {
  kkt_matrix.disp(os);
  return os;
}

} // namespace robotoc 