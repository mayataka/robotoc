#ifndef IDOCP_TERMINAL_OCP_HXX_
#define IDOCP_TERMINAL_OCP_HXX_

#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/ocp/riccati_factorizer.hpp"

#include <cassert>

namespace idocp {

template <typename VectorType>
inline void TerminalOCP::computePrimalDirection(
    const RiccatiFactorization& riccati, 
    const Eigen::MatrixBase<VectorType>& dx0, SplitDirection& d) const {
  static constexpr bool exist_state_constraint = false;
  RiccatiFactorizer::computeStateDirection(riccati, dx0, d, 
                                           exist_state_constraint);
  RiccatiFactorizer::computeCostateDirection(riccati, d, 
                                             exist_state_constraint);
}

} // namespace idocp

#endif // IDOCP_TERMINAL_OCP_HXX_ 