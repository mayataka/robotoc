#ifndef IDOCP_KKT_MATRIX_HPP_
#define IDOCP_KKT_MATRIX_HPP_

#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/ocp/split_switching_constraint_jacobian.hpp"

namespace idocp {

///
/// @typedef KKTMatrix 
/// @brief The KKT matrix of the (hybrid) optimal control problem. 
///
using KKTMatrix = hybrid_container<SplitKKTMatrix, ImpulseSplitKKTMatrix, 
                                   SplitSwitchingConstraintJacobian>;

} // namespace idocp

#endif // IDOCP_KKT_MATRIX_HPP_ 