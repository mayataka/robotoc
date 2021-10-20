#include "idocp/ocp/split_solution.hpp"


namespace idocp {

void SplitSolution::disp(std::ostream& os) const {
  os << "split solution:" << std::endl;
  os << "  q = " << q.transpose() << std::endl;
  os << "  v = " << v.transpose() << std::endl;
  os << "  u = " << u.transpose() << std::endl;
  os << "  a = " << a.transpose() << std::endl;
  if (dimf_ > 0) {
    os << "  f = " << f_stack().transpose() << std::endl;
  }
  os << "  lmd = " << lmd.transpose() << std::endl;
  os << "  gmm = " << gmm.transpose() << std::endl;
  os << "  beta = " << beta.transpose() << std::endl;
  if (dimf_ > 0) {
    os << "  mu = " << mu_stack().transpose() << std::endl;
  }
  if (dimi_ > 0) {
    os << "  xi = " << xi_stack().transpose() << std::flush;
  }
}


std::ostream& operator<<(std::ostream& os, const SplitSolution& s) {
  s.disp(os);
  return os;
}

} // namespace idocp 