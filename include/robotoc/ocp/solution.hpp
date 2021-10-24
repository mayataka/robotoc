#ifndef ROBOTOC_SOLUTION_HPP_
#define ROBOTOC_SOLUTION_HPP_

#include <iostream>

#include "robotoc/hybrid/hybrid_container.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/impulse/impulse_split_solution.hpp"


namespace robotoc {

///
/// @typedef Solution
/// @brief Solution to the (hybrid) optimal control problem. 
///
using Solution = hybrid_container<SplitSolution, ImpulseSplitSolution>;

std::ostream& operator<<(std::ostream& os, const Solution& s);

} // namespace robotoc

#endif // ROBOTOC_SOLUTION_HPP_ 