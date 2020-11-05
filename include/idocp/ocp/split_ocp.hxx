#ifndef IDOCP_SPLIT_OCP_HXX_
#define IDOCP_SPLIT_OCP_HXX_

#include "idocp/ocp/split_ocp.hpp"

#include <cassert>

namespace idocp {

template <typename SplitDirectionType>
inline void SplitOCP::computeCondensedDualDirection(
    Robot& robot, const double dtau, const SplitDirectionType& d_next, 
    SplitDirection& d) {
  assert(dtau > 0);
  contact_dynamics_.computeCondensedDualDirection(robot, dtau, kkt_matrix_,
                                                  kkt_residual_, 
                                                  d_next.dgmm(), d);

}

} // namespace idocp

#endif // IDOCP_SPLIT_OCP_HXX_