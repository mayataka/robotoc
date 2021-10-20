#ifndef IDOCP_DIRECTION_HPP_
#define IDOCP_DIRECTION_HPP_

#include <iostream>

#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"


namespace idocp {

///
/// @typedef Direction
/// @brief Newton direction of the solution to the (hybrid) optimal control 
/// problem. 
///
using Direction = hybrid_container<SplitDirection, ImpulseSplitDirection>;

std::ostream& operator<<(std::ostream& os, const Direction& d);

} // namespace idocp

#endif // IDOCP_DIRECTION_HPP_ 