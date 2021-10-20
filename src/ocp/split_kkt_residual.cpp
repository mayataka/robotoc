#include "idocp/ocp/split_kkt_residual.hpp"


namespace idocp {

void SplitKKTResidual::disp(std::ostream& os) const {
  os << "split KKT residual:" << std::endl;
  os << "  Fq = " << Fq().transpose() << std::endl;
  os << "  Fv = " << Fv().transpose() << std::endl;
  os << "  lq = " << lq().transpose() << std::endl;
  os << "  lv = " << lv().transpose() << std::endl;
  os << "  lu = " << lu.transpose() << std::endl;
  os << "  la = " << la.transpose() << std::endl;
  if (dimf_ > 0) {
    os << "  lf = " << lf().transpose() << std::endl;
  }
  os << "  h = " << h << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const SplitKKTResidual& kkt_residual) {
  kkt_residual.disp(os);
  return os;
}

} // namespace idocp 