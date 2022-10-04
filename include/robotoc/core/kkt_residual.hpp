#ifndef ROBOTOC_KKT_RESIDUAL_HPP_
#define ROBOTOC_KKT_RESIDUAL_HPP_

#include <iostream>

#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/utils/aligned_vector.hpp"


namespace robotoc {

///
/// @typedef KKTResidual 
/// @brief The KKT residual of the optimal control problem. 
///
using KKTResidual = aligned_vector<SplitKKTResidual>;

std::ostream& operator<<(std::ostream& os, const KKTResidual& kkt_residual);

} // namespace robotoc

#endif // ROBOTOC_KKT_RESIDUAL_HPP_ 