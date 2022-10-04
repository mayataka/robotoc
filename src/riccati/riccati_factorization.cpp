#include "robotoc/riccati/riccati_factorization.hpp"


namespace robotoc {

std::ostream& operator<<(std::ostream& os, 
                         const RiccatiFactorization& riccati) {
  os << "RiccatiFactorization:" << "\n";
  for (int i=0; i<riccati.size(); ++i) {
    os << "i: " << i << "\n";
    os << riccati[i] << "\n";
  }
  return os;
}

} // namespace robotoc 