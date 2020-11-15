#ifndef IDOCP_SPLIT_IMPULSE_OCP_HXX_
#define IDOCP_SPLIT_IMPULSE_OCP_HXX_

#include "idocp/impulse/split_impulse_ocp.hpp"

#include <cassert>

namespace idocp {

template <typename MatrixType1, typename MatrixType2>
inline void SplitImpulseOCP::backwardStateConstraintFactorization(
    const Eigen::MatrixBase<MatrixType1>& T_next,
    const Eigen::MatrixBase<MatrixType2>& T) const {
  riccati_factorizer_.backwardStateConstraintFactorization(
      kkt_matrix_, T_next, T);
}


template <typename MatrixType, typename VectorType>
inline void SplitImpulseOCP::getStateConstraintFactorization(
    const Eigen::MatrixBase<MatrixType>& T,
    const Eigen::MatrixBase<VectorType>& e) const {
}


template <typename VectorType>
inline void SplitImpulseOCP::computePrimalDirection(
    Robot& robot, const RiccatiFactorization& riccati, 
    const ImpulseSplitSolution& s, const Eigen::MatrixBase<VectorType>& dx0, 
    ImpulseSplitDirection& d) {
  riccati_factorizer_.computeStateDirection(riccati, dx0, d);
  riccati_factorizer_.computeCostateDirection(riccati, d);
  impulse_dynamics_.computeCondensedPrimalDirection(robot, d);
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, s, d);
}

} // namespace idocp

#endif // IDOCP_SPLIT_IMPULSE_OCP_HXX_ 