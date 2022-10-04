#include "robotoc/riccati/split_riccati_factorization.hpp"


namespace robotoc {

void SplitRiccatiFactorization::disp(std::ostream& os) const {
  os << "split Riccati factorization:" << "\n";
  os << "  P = " << "\n" << P << "\n";
  os << "  s = " << s.transpose() << "\n";
  os << "  psi_x = " << psi_x.transpose() << "\n";
  os << "  psi_u = " << psi_u.transpose() << "\n";
  os << "  Psi = " << Psi.transpose() << "\n";
  os << "  phi_x = " << phi_x.transpose() << "\n";
  os << "  phi_u = " << phi_u.transpose() << "\n";
  os << "  Phi = " << Phi.transpose() << "\n";
  os << "  xi = " << xi << "\n";
  os << "  chi = " << chi << "\n";
  os << "  rho = " << rho << "\n";
  os << "  eta = " << eta << "\n";
  os << "  iota = " << iota << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const SplitRiccatiFactorization& riccati) {
  riccati.disp(os);
  return os;
}

} // namespace robotoc 