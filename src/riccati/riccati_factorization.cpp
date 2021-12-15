#include "robotoc/riccati/riccati_factorization.hpp"


namespace robotoc {

std::ostream& operator<<(std::ostream& os, 
                         const RiccatiFactorization& riccati) {
  os << "Riccati factorization:" << std::endl;
  riccati.disp(os);
  return os;
}

} // namespace robotoc 