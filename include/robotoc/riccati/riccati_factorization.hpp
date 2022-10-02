#ifndef ROBOTOC_RICCATI_FACTORIZATION_HPP_
#define ROBOTOC_RICCATI_FACTORIZATION_HPP_

#include "Eigen/Core"

#include "robotoc/riccati/split_riccati_factorization.hpp"
#include "robotoc/utils/aligned_vector.hpp"


namespace robotoc {

///
/// @typedef RiccatiFactorization
/// @brief Riccati factorization matices of the LQR subproblem. 
///
using RiccatiFactorization = aligned_vector<SplitRiccatiFactorization>;

std::ostream& operator<<(std::ostream& os, const RiccatiFactorization& riccati);

} // namespace robotoc 

#endif // ROBOTOC_RICCATI_FACTORIZATION_HPP_ 