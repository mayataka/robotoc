#ifndef ROBOTOC_DIRECTION_HPP_
#define ROBOTOC_DIRECTION_HPP_

#include <iostream>

#include "robotoc/hybrid/hybrid_container.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_direction.hpp"


namespace robotoc {

///
/// @typedef Direction
/// @brief Newton direction of the solution to the (hybrid) optimal control 
/// problem. 
///
using Direction = hybrid_container<SplitDirection, SplitDirection>;

std::ostream& operator<<(std::ostream& os, const Direction& d);

} // namespace robotoc

#endif // ROBOTOC_DIRECTION_HPP_ 