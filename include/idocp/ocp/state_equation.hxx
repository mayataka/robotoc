#ifndef IDOCP_STATE_EQUATION_HXX_
#define IDOCP_STATE_EQUATION_HXX_

namespace idocp {
  
inline StateEquation::StateEquation(const Robot& robot) 
  : dsubtract_dq_() {
  if (robot.has_floating_base()) {
    dsubtract_dq_.resize(robot.dimv(), robot.dimv());
    dsubtract_dq_.setZero();
  }
}


inline StateEquation::StateEquation() 
  : dsubtract_dq_() {
}


inline StateEquation::~StateEquation() {
}


inline void StateEquation::linearizeForwardEuler(
    Robot& robot, const double dtau, const SplitSolution& s, 
    const SplitSolution& s_next, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  robot.subtractConfiguration(s.q, s_next.q, kkt_residual.Fq());
  kkt_residual.Fq().noalias() += dtau * s.v;
  kkt_residual.Fv() = s.v + dtau * s.a - s_next.v;
  kkt_residual.lq().noalias() += s_next.lmd - s.lmd;
  kkt_residual.lv().noalias() += dtau * s_next.lmd + s_next.gmm - s.gmm;
  kkt_residual.la().noalias() += dtau * s_next.gmm;
  robot.dIntegrateConfiguration(s.q, s.v, dtau, kkt_matrix.Fqq, kkt_matrix.Fqv);
}


template <typename ConfigVectorType1, typename TangentVectorType1, 
          typename TangentVectorType2, typename TangentVectorType3, 
          typename ConfigVectorType2>
inline void StateEquation::linearizeBackwardEuler(
    Robot& robot, const double dtau, 
    const Eigen::MatrixBase<ConfigVectorType1>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType1>& v_prev, 
    const SplitSolution& s, 
    const Eigen::MatrixBase<TangentVectorType2>& lmd_next, 
    const Eigen::MatrixBase<TangentVectorType3>& gmm_next, 
    const Eigen::MatrixBase<ConfigVectorType2>& q_next, 
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  assert(lmd_next.size() == robot.dimv());
  assert(gmm_next.size() == robot.dimv());
  assert(q_next.size() == robot.dimq());
  robot.subtractConfiguration(q_prev, s.q, kkt_residual.Fq());
  kkt_residual.Fq().noalias() += dtau * s.v;
  kkt_residual.Fv() = v_prev - s.v + dtau * s.a;
  if (robot.has_floating_base()) {
    robot.dSubtractdConfigurationMinus(q_prev, s.q, kkt_matrix.Fqq);
    robot.dSubtractdConfigurationPlus(s.q, q_next, dsubtract_dq_);
    kkt_residual.lq().noalias() 
        += dsubtract_dq_.transpose() * lmd_next
            + kkt_matrix.Fqq.transpose() * s.lmd;
  }
  else {
    kkt_residual.lq().noalias() += lmd_next - s.lmd;
  }
  kkt_residual.lv().noalias() += dtau * s.lmd - s.gmm + gmm_next;
  kkt_residual.la().noalias() += dtau * s.gmm;
}


template <typename ConfigVectorType, typename TangentVectorType>
inline void StateEquation::linearizeBackwardEulerTerminal(
    Robot& robot, const double dtau, 
    const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  robot.subtractConfiguration(q_prev, s.q, kkt_residual.Fq());
  kkt_residual.Fq().noalias() += dtau * s.v;
  kkt_residual.Fv() = v_prev - s.v + dtau * s.a;
  if (robot.has_floating_base()) {
    robot.dSubtractdConfigurationMinus(q_prev, s.q, kkt_matrix.Fqq);
    kkt_residual.lq().noalias() 
        += kkt_matrix.Fqq.transpose() * s.lmd;
  }
  else {
    kkt_residual.lq().noalias() -= s.lmd;
  }
  kkt_residual.lv().noalias() += dtau * s.lmd - s.gmm;
  kkt_residual.la().noalias() += dtau * s.gmm;
}


inline double StateEquation::violationL1Norm(
    const KKTResidual& kkt_residual) const {
  return kkt_residual.Fx().lpNorm<1>();
}


template <typename ConfigVectorType, typename TangentVectorType1, 
          typename TangentVectorType2, typename TangentVectorType3>
inline double StateEquation::computeForwardEulerViolationL1Norm(
    const Robot& robot, const double step_size, const double dtau, 
    const SplitSolution& s, const Eigen::MatrixBase<ConfigVectorType>& q_next, 
    const Eigen::MatrixBase<TangentVectorType1>& v_next, 
    const Eigen::MatrixBase<TangentVectorType2>& dq_next, 
    const Eigen::MatrixBase<TangentVectorType3>& dv_next, 
    KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  assert(q_next.size() == robot.dimq());
  assert(v_next.size() == robot.dimv());
  assert(dq_next.size() == robot.dimv());
  assert(dv_next.size() == robot.dimv());
  robot.subtractConfiguration(s.q, q_next, kkt_residual.Fq());
  kkt_residual.Fq().noalias() += dtau * s.v - step_size * dq_next;
  kkt_residual.Fv() = s.v + dtau * s.a - v_next - step_size * dv_next;
  return kkt_residual.Fx().lpNorm<1>();
}


template <typename ConfigVectorType, typename TangentVectorType>
inline double StateEquation::computeBackwardEulerViolationL1Norm(
    const Robot& robot, const double dtau, 
    const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType>& v_prev, const SplitSolution& s, 
    KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  robot.subtractConfiguration(q_prev, s.q, kkt_residual.Fq());
  kkt_residual.Fq().noalias() += dtau * s.v;
  kkt_residual.Fv() = v_prev - s.v + dtau * s.a;
  return kkt_residual.Fx().lpNorm<1>();
}


template <typename ConfigVectorType, typename TangentVectorType1, 
          typename TangentVectorType2, typename TangentVectorType3>
inline double StateEquation::computeBackwardEulerViolationL1Norm(
    const Robot& robot, const double step_size, const double dtau, 
    const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const Eigen::MatrixBase<TangentVectorType1>& v_prev, 
    const Eigen::MatrixBase<TangentVectorType2>& dq_prev, 
    const Eigen::MatrixBase<TangentVectorType3>& dv_prev, 
    const SplitSolution& s, KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  assert(dq_prev.size() == robot.dimv());
  assert(dv_prev.size() == robot.dimv());
  robot.subtractConfiguration(q_prev, s.q, kkt_residual.Fq());
  kkt_residual.Fq().noalias() += dtau * s.v + step_size * dq_prev;
  kkt_residual.Fv() = v_prev - s.v + dtau * s.a + step_size * dv_prev;
  return kkt_residual.Fx().lpNorm<1>();
}

} // namespace idocp 

#endif // IDOCP_STATE_EQUATION_HXX_