#include "robotoc/riccati/split_constrained_riccati_factorization.hpp"


namespace robotoc {

void SplitConstrainedRiccatiFactorization::disp(std::ostream& os) const {
}


std::ostream& operator<<(std::ostream& os, 
                         const SplitConstrainedRiccatiFactorization& c_riccati) {
  c_riccati.disp(os);
  return os;
}

} // namespace robotoc