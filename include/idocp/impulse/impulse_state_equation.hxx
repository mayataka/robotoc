#ifndef IDOCP_IMPULSE_STATE_EQUATION_HXX_
#define IDOCP_IMPULSE_STATE_EQUATION_HXX_

#include "idocp/impulse/impulse_state_equation.hpp"

#include <cassert>

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
    kkt_residual.lq().template head<6>().noalias() 
        += kkt_matrix.Fqq().template topLeftCorner<6, 6>().transpose() 
              * s_next.lmd.template head<6>();
    kkt_residual.lq().template head<6>().noalias() 
        += kkt_matrix.Fqq_prev.template topLeftCorner<6, 6>().transpose() 
              * s.lmd.template head<6>();
    kkt_residual.lq().tail(robot.dimv()-6).noalias() 
        += s_next.lmd.tail(robot.dimv()-6) - s.lmd.tail(robot.dimv()-6);
  }
  else {
    kkt_matrix.Fqq().diagonal().fill(1.);
    kkt_residual.lq().noalias() += s_next.lmd - s.lmd;
  }
  kkt_residual.lv().noalias() += s_next.gmm - s.gmm;
  kkt_residual.ldv.noalias() += s_next.gmm;
}


template <typename ConfigVectorType>
inline void condenseImpulseForwardEuler(
    Robot& robot, const ImpulseSplitSolution& s, 
    const Eigen::MatrixBase<ConfigVectorType>& q_next, 
    ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual) {
  if (robot.hasFloatingBase()) {
    robot.dSubtractdConfigurationInverse(kkt_matrix.Fqq_prev, 
                                         kkt_matrix.Fqq_prev_inv);
    robot.dSubtractdConfigurationMinus(s.q, q_next, kkt_matrix.Fqq_prev);
    robot.dSubtractdConfigurationInverse(kkt_matrix.Fqq_prev, 
                                         kkt_matrix.Fqq_inv);
    kkt_matrix.Fqq_prev.template topLeftCorner<6, 6>() 
        = kkt_matrix.Fqq().template topLeftCorner<6, 6>();
    kkt_residual.Fq_tmp.template head<6>() 
        = kkt_residual.Fq().template head<6>();
    kkt_matrix.Fqq().template topLeftCorner<6, 6>().noalias()
        = - kkt_matrix.Fqq_inv 
            * kkt_matrix.Fqq_prev.template topLeftCorner<6, 6>();
    kkt_residual.Fq().template head<6>().noalias()
        = - kkt_matrix.Fqq_inv * kkt_residual.Fq_tmp.template head<6>();
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


inline double l1NormStateEuqationResidual(
    const ImpulseSplitKKTResidual& kkt_residual) {
  assert(kkt_residual.isDimensionConsistent());
  return kkt_residual.Fx.lpNorm<1>();
}


inline double squaredNormStateEuqationResidual(
    const ImpulseSplitKKTResidual& kkt_residual) {
  assert(kkt_residual.isDimensionConsistent());
  return kkt_residual.Fx.squaredNorm();
}

} // namespace stateequation
} // namespace idocp 

#endif // IDOCP_IMPULSE_STATE_EQUATION_HXX_