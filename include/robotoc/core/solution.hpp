#ifndef ROBOTOC_SOLUTION_HPP_
#define ROBOTOC_SOLUTION_HPP_

#include <iostream>

#include "robotoc/core/split_solution.hpp"
#include "robotoc/utils/aligned_vector.hpp"


namespace robotoc {

///
/// @typedef Solution
/// @brief Solution to the optimal control problem. 
///
using Solution = aligned_vector<SplitSolution>;

std::ostream& operator<<(std::ostream& os, const Solution& s);

} // namespace robotoc

#endif // ROBOTOC_SOLUTION_HPP_ 