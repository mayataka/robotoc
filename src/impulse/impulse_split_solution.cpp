#include "idocp/impulse/impulse_split_solution.hpp"


namespace idocp {

void ImpulseSplitSolution::disp(std::ostream& os) const {
  os << "impulse split solution:" << std::endl;
  os << "  q = " << q.transpose() << std::endl;
  os << "  v = " << v.transpose() << std::endl;
  os << "  dv = " << dv.transpose() << std::endl;
  if (dimi_ > 0) {
    os << "  f = " << f_stack().transpose() << std::endl;
  }
  os << "  lmd = " << lmd.transpose() << std::endl;
  os << "  gmm = " << gmm.transpose() << std::endl;
  os << "  beta = " << beta.transpose() << std::endl;
  if (dimi_ > 0) {
    os << "  mu = " << mu_stack().transpose() << std::endl;
  }
}


std::ostream& operator<<(std::ostream& os, const ImpulseSplitSolution& s) {
  s.disp(os);
  return os;
}

} // namespace idocp 