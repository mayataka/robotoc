#ifndef ROBOTOC_SOLUTION_HPP_
#define ROBOTOC_SOLUTION_HPP_

#include <iostream>

#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/hybrid_container.hpp"


namespace robotoc {

///
/// @typedef Solution
/// @brief Solution to the (hybrid) optimal control problem. 
///
using Solution = hybrid_container<SplitSolution, SplitSolution>;

std::ostream& operator<<(std::ostream& os, const Solution& s);

} // namespace robotoc

#endif // ROBOTOC_SOLUTION_HPP_ 