#include "robotoc/impulse/impulse_split_kkt_residual.hpp"


namespace robotoc {

void ImpulseSplitKKTResidual::disp(std::ostream& os) const {
  os << "impulse split KKT residual:" << std::endl;
  os << "  Fq = " << Fq().transpose() << std::endl;
  os << "  Fv = " << Fv().transpose() << std::endl;
  os << "  lq = " << lq().transpose() << std::endl;
  os << "  lv = " << lv().transpose() << std::endl;
  os << "  ldv = " << ldv.transpose() << std::endl;
  if (dimi_ > 0) {
    os << "  lf = " << lf().transpose() << std::endl;
  }
  os << "  kkt_error = " << kkt_error << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const ImpulseSplitKKTResidual& kkt_residual) {
  kkt_residual.disp(os);
  return os;
}

} // namespace robotoc 