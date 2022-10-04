#ifndef ROBOTOC_DIRECTION_HPP_
#define ROBOTOC_DIRECTION_HPP_

#include <iostream>

#include "robotoc/core/split_direction.hpp"
#include "robotoc/utils/aligned_vector.hpp"


namespace robotoc {

///
/// @typedef Direction
/// @brief Newton direction of the solution to the optimal control problem. 
///
using Direction = aligned_vector<SplitDirection>;

std::ostream& operator<<(std::ostream& os, const Direction& d);

} // namespace robotoc

#endif // ROBOTOC_DIRECTION_HPP_ 