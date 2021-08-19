#ifndef IDOCP_SOLUTION_HPP_
#define IDOCP_SOLUTION_HPP_

#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"

namespace idocp {

///
/// @typedef Solution
/// @brief Solution to the (hybrid) optimal control problem. 
///
using Solution = hybrid_container<SplitSolution, ImpulseSplitSolution>;

} // namespace idocp

#endif // IDOCP_SOLUTION_HPP_ 