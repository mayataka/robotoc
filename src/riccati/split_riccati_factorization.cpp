#include "robotoc/riccati/split_riccati_factorization.hpp"


namespace robotoc {

void SplitRiccatiFactorization::disp(std::ostream& os) const {
  os << "split Riccati factorization:" << std::endl;
  os << "  P = " << P << std::endl;
  os << "  s = " << s.transpose() << std::endl;
  os << "  psi_x = " << psi_x.transpose() << std::endl;
  os << "  psi_u = " << psi_u.transpose() << std::endl;
  os << "  Psi = " << Psi.transpose() << std::endl;
  os << "  phi_x = " << phi_x.transpose() << std::endl;
  os << "  phi_u = " << phi_u.transpose() << std::endl;
  os << "  Phi = " << Phi.transpose() << std::endl;
  os << "  xi = " << xi << std::endl;
  os << "  chi = " << chi << std::endl;
  os << "  rho = " << rho << std::endl;
  os << "  eta = " << eta << std::endl;
  os << "  iota = " << iota << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const SplitRiccatiFactorization& riccati) {
  riccati.disp(os);
  return os;
}

} // namespace robotoc 