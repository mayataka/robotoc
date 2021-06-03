#ifndef IDOCP_UNCONSTR_STATE_EQUATION_HXX_
#define IDOCP_UNCONSTR_STATE_EQUATION_HXX_

#include "idocp/unconstr/unconstr_state_equation.hpp"

#include <cassert>

namespace idocp {
namespace unconstr {
namespace stateequation {

inline void linearizeForwardEuler(const double dt, const SplitSolution& s, 
                                  const SplitSolution& s_next, 
                                  SplitKKTMatrix& kkt_matrix, 
                                  SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  computeForwardEulerResidual(dt, s, s_next.q, s_next.v, kkt_residual);
  kkt_residual.lq().noalias() += s_next.lmd - s.lmd;
  kkt_residual.lv().noalias() += dt * s_next.lmd + s_next.gmm - s.gmm;
  kkt_residual.la.noalias()   += dt * s_next.gmm;
}


template <typename ConfigVectorType, typename TangentVectorType>
inline void linearizeBackwardEuler(
    const double dt, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const SplitSolution& s, const SplitSolution& s_next,
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  computeBackwardEulerResidual(dt, q_prev, v_prev, s, kkt_residual);
  kkt_residual.lq().noalias() += s_next.lmd - s.lmd;
  kkt_residual.lv().noalias() += dt * s.lmd - s.gmm + s_next.gmm;
  kkt_residual.la.noalias()   += dt * s.gmm;
}


template <typename ConfigVectorType, typename TangentVectorType>
inline void linearizeBackwardEulerTerminal(
    const double dt, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  computeBackwardEulerResidual(dt, q_prev, v_prev, s, kkt_residual);
  kkt_residual.lq().noalias() -= s.lmd;
  kkt_residual.lv().noalias() += dt * s.lmd - s.gmm;
  kkt_residual.la.noalias()   += dt * s.gmm;
}


template <typename ConfigVectorType, typename TangentVectorType>
inline void computeForwardEulerResidual(
    const double dt, const SplitSolution& s, 
    const Eigen::MatrixBase<ConfigVectorType>& q_next, 
    const Eigen::MatrixBase<TangentVectorType>& v_next, 
    SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  kkt_residual.Fq() = s.q + dt * s.v - q_next;
  kkt_residual.Fv() = s.v + dt * s.a - v_next;
}


template <typename ConfigVectorType, typename TangentVectorType>
inline void computeBackwardEulerResidual(
    const double dt, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  kkt_residual.Fq() = q_prev - s.q + dt * s.v;
  kkt_residual.Fv() = v_prev - s.v + dt * s.a;
}


inline double l1NormStateEuqationResidual(
    const SplitKKTResidual& kkt_residual) {
  return kkt_residual.Fx.lpNorm<1>();
}


inline double squaredNormStateEuqationResidual(
    const SplitKKTResidual& kkt_residual) {
  return kkt_residual.Fx.squaredNorm();
}

} // namespace stateequation 
} // namespace unconstr 
} // namespace idocp 

#endif // IDOCP_UNCONSTR_STATE_EQUATION_HXX_ 