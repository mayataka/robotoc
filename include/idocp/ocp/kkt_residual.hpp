#ifndef IDOCP_KKT_RESIDUAL_HPP_
#define IDOCP_KKT_RESIDUAL_HPP_

#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/ocp/split_switching_constraint_residual.hpp"

namespace idocp {

///
/// @typedef KKTResidual 
/// @brief The KKT residual of the (hybrid) optimal control problem. 
///
using KKTResidual = hybrid_container<SplitKKTResidual, ImpulseSplitKKTResidual, 
                                     SplitSwitchingConstraintResidual>;

} // namespace idocp

#endif // IDOCP_KKT_RESIDUAL_HPP_ 