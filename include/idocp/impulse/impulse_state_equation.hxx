#ifndef IDOCP_IMPULSE_STATE_EQUATION_HXX_
#define IDOCP_IMPULSE_STATE_EQUATION_HXX_

#include "idocp/impulse/impulse_state_equation.hpp"

#include <assert.h>

namespace idocp {
namespace stateequation {

template <typename ConfigVectorType>
inline void LinearizeImpulseForwardEuler(
    const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const ImpulseSplitSolution& s, const SplitSolution& s_next, 
    ImpulseKKTMatrix& kkt_matrix, ImpulseKKTResidual& kkt_residual) {
  ComputeImpulseForwardEulerResidual(robot, s, s_next.q, s_next.v, kkt_residual);
  if (robot.has_floating_base()) {
    robot.dSubtractdConfigurationPlus(s.q, s_next.q, kkt_matrix.Fqq);
    robot.dSubtractdConfigurationMinus(q_prev, s.q, kkt_matrix.Fqq_prev);
    kkt_residual.lq().noalias() += kkt_matrix.Fqq.transpose() * s_next.lmd 
                                    + kkt_matrix.Fqq_prev.transpose() * s.lmd;
  }
  else {
    kkt_residual.lq().noalias() += s_next.lmd - s.lmd;
  }
  kkt_residual.lv().noalias() += s_next.gmm - s.gmm;
  kkt_residual.ldv.noalias() += s_next.gmm;
}


template <typename ConfigVectorType, typename TangentVectorType>
inline void LinearizeImpulseBackwardEuler(
    const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const ImpulseSplitSolution& s, const SplitSolution& s_next,
    ImpulseKKTMatrix& kkt_matrix, ImpulseKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  ComputeImpulseBackwardEulerResidual(robot, q_prev, v_prev, s, kkt_residual);
  if (robot.has_floating_base()) {
    robot.dSubtractdConfigurationMinus(q_prev, s.q, kkt_matrix.Fqq);
    robot.dSubtractdConfigurationPlus(s.q, s_next.q, kkt_matrix.Fqq_prev);
    kkt_residual.lq().noalias() 
        += kkt_matrix.Fqq_prev.transpose() * s_next.lmd
            + kkt_matrix.Fqq.transpose() * s.lmd;
  }
  else {
    kkt_residual.lq().noalias() += s_next.lmd - s.lmd;
  }
  kkt_residual.lv().noalias() += - s.gmm + s_next.gmm;
  kkt_residual.ldv.noalias() += s.gmm;
}


template <typename ConfigVectorType, typename TangentVectorType>
inline void LinearizeImpulseBackwardEulerTerminal(
    const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const ImpulseSplitSolution& s, ImpulseKKTMatrix& kkt_matrix, 
    ImpulseKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  ComputeImpulseBackwardEulerResidual(robot, q_prev, v_prev, s, kkt_residual);
  if (robot.has_floating_base()) {
    robot.dSubtractdConfigurationMinus(q_prev, s.q, kkt_matrix.Fqq);
    kkt_residual.lq().noalias() 
        += kkt_matrix.Fqq.transpose() * s.lmd;
  }
  else {
    kkt_residual.lq().noalias() -= s.lmd;
  }
  kkt_residual.lv().noalias() -= s.gmm;
  kkt_residual.ldv.noalias() += s.gmm;
}


template <typename ConfigVectorType, typename TangentVectorType>
inline void ComputeImpulseForwardEulerResidual(
    const Robot& robot, const ImpulseSplitSolution& s, 
    const Eigen::MatrixBase<ConfigVectorType>& q_next, 
    const Eigen::MatrixBase<TangentVectorType>& v_next, 
    ImpulseKKTResidual& kkt_residual) {
  assert(q_next.size() == robot.dimq());
  assert(v_next.size() == robot.dimv());
  robot.subtractConfiguration(s.q, q_next, kkt_residual.Fq());
  kkt_residual.Fv() = s.v + s.dv - v_next;
}


template <typename ConfigVectorType, typename TangentVectorType1, 
          typename TangentVectorType2, typename TangentVectorType3>
inline void ComputeImpulseForwardEulerResidual(
    const Robot& robot, const double step_size, const ImpulseSplitSolution& s, 
    const Eigen::MatrixBase<ConfigVectorType>& q_next, 
    const Eigen::MatrixBase<TangentVectorType1>& v_next, 
    const Eigen::MatrixBase<TangentVectorType2>& dq_next, 
    const Eigen::MatrixBase<TangentVectorType3>& dv_next, 
    ImpulseKKTResidual& kkt_residual) {
  assert(q_next.size() == robot.dimq());
  assert(v_next.size() == robot.dimv());
  assert(dq_next.size() == robot.dimv());
  assert(dv_next.size() == robot.dimv());
  robot.subtractConfiguration(s.q, q_next, kkt_residual.Fq());
  kkt_residual.Fq().noalias() -= step_size * dq_next;
  kkt_residual.Fv() = s.v + s.dv - v_next - step_size * dv_next;
}


template <typename ConfigVectorType, typename TangentVectorType>
inline void ComputeImpulseBackwardEulerResidual(
    const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const ImpulseSplitSolution& s, ImpulseKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  robot.subtractConfiguration(q_prev, s.q, kkt_residual.Fq());
  kkt_residual.Fv() = v_prev - s.v + s.dv;
}


template <typename ConfigVectorType, typename TangentVectorType1, 
          typename TangentVectorType2, typename TangentVectorType3>
inline void ComputeImpulseBackwardEulerResidual(
    const Robot& robot, const double step_size, 
    const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType1>& v_prev, 
    const Eigen::MatrixBase<TangentVectorType2>& dq_prev, 
    const Eigen::MatrixBase<TangentVectorType3>& dv_prev, 
    const ImpulseSplitSolution& s, ImpulseKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  assert(dq_prev.size() == robot.dimv());
  assert(dv_prev.size() == robot.dimv());
  robot.subtractConfiguration(q_prev, s.q, kkt_residual.Fq());
  kkt_residual.Fq().noalias() += step_size * dq_prev;
  kkt_residual.Fv() = v_prev - s.v + s.dv + step_size * dv_prev;
}


inline double L1NormStateEuqationResidual(
    const ImpulseKKTResidual& kkt_residual) {
  return kkt_residual.Fx().lpNorm<1>();
}


inline double SquaredNormStateEuqationResidual(
    const ImpulseKKTResidual& kkt_residual) {
  return kkt_residual.Fx().squaredNorm();
}

} // namespace stateequation
} // namespace idocp 

#endif // IDOCP_IMPULSE_STATE_EQUATION_HXX_