#ifndef IDOCP_IMPULSE_STATE_EQUATION_HXX_
#define IDOCP_IMPULSE_STATE_EQUATION_HXX_

#include "idocp/impulse/impulse_state_equation.hpp"

#include <assert.h>

namespace idocp {
namespace stateequation {

template <typename ConfigVectorType>
inline void linearizeImpulseForwardEuler(
    const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const ImpulseSplitSolution& s, const SplitSolution& s_next, 
    ImpulseSplitKKTMatrix& kkt_matrix, ImpulseSplitKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  computeImpulseForwardEulerResidual(robot, s, s_next.q, s_next.v, kkt_residual);
  if (robot.hasFloatingBase()) {
    robot.dSubtractdConfigurationPlus(s.q, s_next.q, kkt_matrix.Fqq());
    robot.dSubtractdConfigurationMinus(q_prev, s.q, kkt_matrix.Fqq_prev);
    kkt_residual.lq().noalias() += kkt_matrix.Fqq().transpose() * s_next.lmd 
                                    + kkt_matrix.Fqq_prev.transpose() * s.lmd;
  }
  else {
    kkt_residual.lq().noalias() += s_next.lmd - s.lmd;
  }
  kkt_residual.lv().noalias() += s_next.gmm - s.gmm;
  kkt_residual.ldv.noalias() += s_next.gmm;
}


inline void condenseImpulseForwardEuler(Robot& robot, 
                                        ImpulseSplitKKTMatrix& kkt_matrix, 
                                        ImpulseSplitKKTResidual& kkt_residual) {
  if (robot.hasFloatingBase()) {
    robot.dSubtractdConfigurationInverse(kkt_matrix.Fqq_prev, 
                                         kkt_matrix.Fqq_prev_inv);
    kkt_matrix.Fqq_prev.template topLeftCorner<6, 6>() 
        = kkt_matrix.Fqq().template topLeftCorner<6, 6>();
    kkt_residual.Fq_prev = kkt_residual.Fq().template head<6>();
    kkt_matrix.Fqq().template topLeftCorner<6, 6>().noalias()
        = - kkt_matrix.Fqq_prev_inv 
            * kkt_matrix.Fqq_prev.template topLeftCorner<6, 6>();
    kkt_residual.Fq().template head<6>().noalias()
        = - kkt_matrix.Fqq_prev_inv * kkt_residual.Fq_prev;
  }
}


template <typename ConfigVectorType, typename TangentVectorType>
inline void linearizeImpulseBackwardEuler(
    const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const ImpulseSplitSolution& s, const SplitSolution& s_next,
    ImpulseSplitKKTMatrix& kkt_matrix, ImpulseSplitKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  computeImpulseBackwardEulerResidual(robot, q_prev, v_prev, s, kkt_residual);
  if (robot.hasFloatingBase()) {
    robot.dSubtractdConfigurationMinus(q_prev, s.q, kkt_matrix.Fqq());
    robot.dSubtractdConfigurationPlus(s.q, s_next.q, kkt_matrix.Fqq_prev);
    kkt_residual.lq().noalias() 
        += kkt_matrix.Fqq_prev.transpose() * s_next.lmd
            + kkt_matrix.Fqq().transpose() * s.lmd;
  }
  else {
    kkt_residual.lq().noalias() += s_next.lmd - s.lmd;
  }
  kkt_residual.lv().noalias() += - s.gmm + s_next.gmm;
  kkt_residual.ldv.noalias() += s.gmm;
}


template <typename ConfigVectorType>
inline void condenseImpulseBackwardEuler(
    Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual) {
  if (robot.hasFloatingBase()) {
    robot.dSubtractdConfigurationPlus(q_prev, s.q, kkt_matrix.Fqq_prev);
    robot.dSubtractdConfigurationInverse(kkt_matrix.Fqq_prev, 
                                         kkt_matrix.Fqq_prev_inv);
    kkt_matrix.Fqq_prev.template topLeftCorner<6, 6>() 
        = kkt_matrix.Fqq().template topLeftCorner<6, 6>();
    kkt_residual.Fq_prev = kkt_residual.Fq().template head<6>();
    kkt_matrix.Fqq().template topLeftCorner<6, 6>().noalias()
        = kkt_matrix.Fqq_prev_inv 
            * kkt_matrix.Fqq_prev.template topLeftCorner<6, 6>();
    kkt_residual.Fq().template head<6>().noalias()
        = kkt_matrix.Fqq_prev_inv * kkt_residual.Fq_prev;
  }
}


template <typename ConfigVectorType, typename TangentVectorType>
inline void computeImpulseForwardEulerResidual(
    const Robot& robot, const ImpulseSplitSolution& s, 
    const Eigen::MatrixBase<ConfigVectorType>& q_next, 
    const Eigen::MatrixBase<TangentVectorType>& v_next, 
    ImpulseSplitKKTResidual& kkt_residual) {
  assert(q_next.size() == robot.dimq());
  assert(v_next.size() == robot.dimv());
  robot.subtractConfiguration(s.q, q_next, kkt_residual.Fq());
  kkt_residual.Fv() = s.v + s.dv - v_next;
}


template <typename ConfigVectorType, typename TangentVectorType>
inline void computeImpulseBackwardEulerResidual(
    const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  robot.subtractConfiguration(q_prev, s.q, kkt_residual.Fq());
  kkt_residual.Fv() = v_prev - s.v + s.dv;
}


inline double l1NormStateEuqationResidual(
    const ImpulseSplitKKTResidual& kkt_residual) {
  return kkt_residual.Fx().lpNorm<1>();
}


inline double squaredNormStateEuqationResidual(
    const ImpulseSplitKKTResidual& kkt_residual) {
  return kkt_residual.Fx().squaredNorm();
}

} // namespace stateequation
} // namespace idocp 

#endif // IDOCP_IMPULSE_STATE_EQUATION_HXX_