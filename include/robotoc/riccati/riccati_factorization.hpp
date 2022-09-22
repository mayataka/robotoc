#ifndef ROBOTOC_RICCATI_FACTORIZATION_HPP_
#define ROBOTOC_RICCATI_FACTORIZATION_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/ocp/hybrid_container.hpp"
#include "robotoc/riccati/split_riccati_factorization.hpp"
#include "robotoc/riccati/split_constrained_riccati_factorization.hpp"


namespace robotoc {

///
/// @typedef RiccatiFactorization
/// @brief Riccati factorization matices of the LQR subproblem. 
///
using RiccatiFactorization = hybrid_container<SplitRiccatiFactorization, 
                                              SplitRiccatiFactorization, 
                                              SplitConstrainedRiccatiFactorization>;

std::ostream& operator<<(std::ostream& os, const RiccatiFactorization& riccati);

} // namespace robotoc 

#endif // ROBOTOC_RICCATI_FACTORIZATION_HPP_ 