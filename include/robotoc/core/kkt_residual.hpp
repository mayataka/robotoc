#ifndef ROBOTOC_KKT_RESIDUAL_HPP_
#define ROBOTOC_KKT_RESIDUAL_HPP_

#include <iostream>

#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/switching_constraint_residual.hpp"
#include "robotoc/core/hybrid_container.hpp"


namespace robotoc {

///
/// @typedef KKTResidual 
/// @brief The KKT residual of the (hybrid) optimal control problem. 
///
using KKTResidual = hybrid_container<SplitKKTResidual, SplitKKTResidual, 
                                     SwitchingConstraintResidual>;

std::ostream& operator<<(std::ostream& os, const KKTResidual& kkt_residual);

} // namespace robotoc

#endif // ROBOTOC_KKT_RESIDUAL_HPP_ 