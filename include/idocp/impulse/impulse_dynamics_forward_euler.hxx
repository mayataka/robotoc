#ifndef IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HXX_
#define IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HXX_

#include "idocp/impulse/impulse_dynamics_forward_euler.hpp"

#include <cassert>

namespace idocp {

inline ImpulseDynamicsForwardEuler::ImpulseDynamicsForwardEuler(
    const Robot& robot) 
  : data_(robot),
    penalty_(0) {
}


inline ImpulseDynamicsForwardEuler::ImpulseDynamicsForwardEuler() 
  : data_(),
    penalty_(0) {
}


inline ImpulseDynamicsForwardEuler::~ImpulseDynamicsForwardEuler() {
}


inline void ImpulseDynamicsForwardEuler::linearizeImpulseDynamics(
    Robot& robot, const ImpulseStatus& impulse_status,  
    const ImpulseSplitSolution& s, ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual, 
    const bool is_state_constraint_valid) {
  setImpulseStatus(impulse_status);
  linearizeInverseImpulseDynamics(robot, impulse_status, s, data_);
  linearizeImpulseVelocityConstraint(robot, impulse_status, data_);
  // augment inverse impulse dynamics constraint
  kkt_residual.lq().noalias() += data_.dImDdq().transpose() * s.beta;
  kkt_residual.ldv.noalias() += data_.dImDddv.transpose() * s.beta;
  // We use an equivalence dmIDdf_().transpose() = - dCdv_() = - dCddv, to avoid
  // redundant calculation of dImDdf_().
  kkt_residual.lf().noalias() -= data_.dCdv() * s.beta;
  // augment impulse velocity constraint
  kkt_residual.lq().noalias() += data_.dCdq().transpose() * s.mu_stack();
  kkt_residual.lv().noalias() += data_.dCdv().transpose() * s.mu_stack();
  // We use an equivalence dCdv_() = dCddv, to avoid redundant calculation.
  kkt_residual.ldv.noalias() += data_.dCdv().transpose() * s.mu_stack();
  // augment impulse position constraint
  if (is_state_constraint_valid) {
    linearizeImpulsePositionConstraint(robot, impulse_status, kkt_matrix, 
                                       kkt_residual);
    kkt_residual.lq().noalias() 
        += s.penalty * kkt_matrix.Pq().transpose() * kkt_residual.P();
    kkt_residual.lq().noalias() 
        += kkt_matrix.Pq().transpose() * s.xi_stack();
    penalty_ = s.penalty;
  }
}


inline void ImpulseDynamicsForwardEuler::linearizeInverseImpulseDynamics(
    Robot& robot, const ImpulseStatus& impulse_status, 
    const ImpulseSplitSolution& s, ImpulseDynamicsForwardEulerData& data) {
  robot.setImpulseForces(impulse_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, data.ImD());
  robot.RNEAImpulseDerivatives(s.q, s.dv, data.dImDdq(), data.dImDddv);
}


inline void ImpulseDynamicsForwardEuler::linearizeImpulseVelocityConstraint(
    Robot& robot, const ImpulseStatus& impulse_status, 
    ImpulseDynamicsForwardEulerData& data) {
  robot.computeImpulseVelocityResidual(impulse_status, data.C());
  robot.computeImpulseVelocityDerivatives(impulse_status, data.dCdq(), 
                                          data.dCdv());
}


inline void ImpulseDynamicsForwardEuler::linearizeImpulsePositionConstraint(
      Robot& robot, const ImpulseStatus& impulse_status, 
      ImpulseSplitKKTMatrix& kkt_matrix, ImpulseSplitKKTResidual& kkt_residual) {
  robot.computeImpulseConditionResidual(impulse_status,
                                        impulse_status.contactPoints(),
                                        kkt_residual.P());
  robot.computeImpulseConditionDerivative(impulse_status, kkt_matrix.Pq());
}


inline void ImpulseDynamicsForwardEuler::condenseImpulseDynamics(
    Robot& robot, const ImpulseStatus& impulse_status, 
    ImpulseSplitKKTMatrix& kkt_matrix, ImpulseSplitKKTResidual& kkt_residual) {
  robot.computeMJtJinv(data_.dImDddv, data_.dCdv(), data_.MJtJinv());
  condensing(robot, impulse_status, data_, kkt_matrix, kkt_residual);
  kkt_matrix.Qqq().noalias() 
      += penalty_ * kkt_matrix.Pq().transpose() * kkt_matrix.Pq();
}


inline void ImpulseDynamicsForwardEuler::condensing(
    const Robot& robot, const ImpulseStatus& impulse_status, 
    ImpulseDynamicsForwardEulerData& data, ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual) {
  const int dimv = robot.dimv();
  const int dimf = impulse_status.dimf();
  data.MJtJinv_dImDCdqv().leftCols(dimv).noalias() 
      = data.MJtJinv() * data.dImDCdq();
  data.MJtJinv_dImDCdqv().topRightCorner(dimv, dimv).noalias() 
      = data.MJtJinv().topRightCorner(dimv, dimf) * data.dCdv();
  data.MJtJinv_dImDCdqv().bottomRightCorner(dimf, dimv).noalias() 
      = data.MJtJinv().bottomRightCorner(dimf, dimf) * data.dCdv();
  data.MJtJinv_ImDC().noalias() = data.MJtJinv() * data.ImDC();
  data.Qdvfqv().topRows(dimv).noalias() 
      = (- kkt_matrix.Qdvdv().diagonal()).asDiagonal() 
          * data.MJtJinv_dImDCdqv().topRows(dimv);
  data.Qdvfqv().bottomRows(dimf).noalias() 
      = - kkt_matrix.Qff() * data.MJtJinv_dImDCdqv().bottomRows(dimf);
  data.ldv() = kkt_residual.ldv;
  data.lf()  = - kkt_residual.lf();
  data.ldv().noalias() 
      -= kkt_matrix.Qdvdv().diagonal().asDiagonal() 
          * data.MJtJinv_ImDC().head(dimv);
  data.lf().noalias() -= kkt_matrix.Qff() * data.MJtJinv_ImDC().tail(dimf);
  kkt_matrix.Qxx().noalias() 
      -= data.MJtJinv_dImDCdqv().transpose() * data.Qdvfqv();
  kkt_residual.lx().noalias() 
      -= data.MJtJinv_dImDCdqv().transpose() * data.ldvf();
  kkt_matrix.Fvq() = - data.MJtJinv_dImDCdqv().topLeftCorner(dimv, dimv);
  kkt_matrix.Fvv() = Eigen::MatrixXd::Identity(dimv, dimv) 
                    - data.MJtJinv_dImDCdqv().topRightCorner(dimv, dimv);
  kkt_residual.Fv().noalias() -= data.MJtJinv_ImDC().head(dimv);
}


inline void ImpulseDynamicsForwardEuler::computeCondensedPrimalDirection(
    const Robot& robot, ImpulseSplitDirection& d) {
  expansionPrimal(robot, data_, d);
}


template <typename VectorType>
inline void ImpulseDynamicsForwardEuler::computeCondensedDualDirection(
    const Robot& robot, const ImpulseSplitKKTMatrix& kkt_matrix, 
    const ImpulseSplitKKTResidual& kkt_residual, 
    const Eigen::MatrixBase<VectorType>& dgmm, ImpulseSplitDirection& d) {
  assert(dgmm.size() == robot.dimv());
  expansionDual(robot, data_, kkt_matrix, kkt_residual, dgmm, d);
}


inline void ImpulseDynamicsForwardEuler::expansionPrimal(
    const Robot& robot, const ImpulseDynamicsForwardEulerData& data, 
    ImpulseSplitDirection& d) {
  d.ddvf().noalias() = - data.MJtJinv_dImDCdqv() * d.dx();
  d.ddvf().noalias() -= data.MJtJinv_ImDC();
  d.df().array() *= -1;
}


template <typename VectorType>
inline void ImpulseDynamicsForwardEuler::expansionDual(
    const Robot& robot, ImpulseDynamicsForwardEulerData& data, 
    const ImpulseSplitKKTMatrix& kkt_matrix, 
    const ImpulseSplitKKTResidual& kkt_residual, 
    const Eigen::MatrixBase<VectorType>& dgmm, ImpulseSplitDirection& d) {
  assert(dgmm.size() == robot.dimv());
  data.ldvf().noalias() += data.Qdvfqv() * d.dx();
  data.ldv().noalias() += dgmm;
  d.dbetamu().noalias() = - data.MJtJinv() * data.ldvf();
}


inline void ImpulseDynamicsForwardEuler::computeImpulseDynamicsResidual(
    Robot& robot, const ImpulseStatus& impulse_status,
    const ImpulseSplitSolution& s, ImpulseSplitKKTResidual& kkt_residual,
    const bool is_state_constraint_valid) {
  setImpulseStatus(impulse_status);
  robot.setImpulseForces(impulse_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, data_.ImD());
  robot.computeImpulseVelocityResidual(impulse_status, data_.C());
  if (is_state_constraint_valid) {
    robot.computeImpulseConditionResidual(impulse_status, 
                                          impulse_status.contactPoints(), 
                                          kkt_residual.P());
  }
}


inline double ImpulseDynamicsForwardEuler::l1NormImpulseDynamicsResidual(
    const ImpulseSplitKKTResidual& kkt_residual,
    const bool is_state_constraint_valid) const {
  return data_.ImDC().lpNorm<1>();
}


inline double ImpulseDynamicsForwardEuler::squaredNormImpulseDynamicsResidual(
    const ImpulseSplitKKTResidual& kkt_residual,
    const bool is_state_constraint_valid) const {
  return (data_.ImDC().squaredNorm() + kkt_residual.P().squaredNorm());
}


inline void ImpulseDynamicsForwardEuler::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  data_.setImpulseStatus(impulse_status);
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HXX_ 