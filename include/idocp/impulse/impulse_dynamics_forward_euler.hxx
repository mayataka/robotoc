#ifndef IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HXX_
#define IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HXX_

#include "idocp/impulse/impulse_dynamics_forward_euler.hpp"

#include <cassert>

namespace idocp {

inline ImpulseDynamicsForwardEuler::ImpulseDynamicsForwardEuler(
    const Robot& robot) 
  : data_(robot) {
}


inline ImpulseDynamicsForwardEuler::ImpulseDynamicsForwardEuler() 
  : data_() {
}


inline ImpulseDynamicsForwardEuler::~ImpulseDynamicsForwardEuler() {
}


inline void ImpulseDynamicsForwardEuler::linearizeImpulseDynamics(
    Robot& robot, const ImpulseStatus& impulse_status,  
    const ImpulseSplitSolution& s, ImpulseKKTMatrix& kkt_matrix, 
    ImpulseKKTResidual& kkt_residual) { 
  setImpulseStatus(impulse_status);
  linearizeInverseImpulseDynamics(robot, impulse_status, s, data_);
  linearizeImpulseVelocityConstraint(robot, impulse_status, data_);
  linearizeImpulsePositionConstraint(robot, impulse_status, kkt_matrix, 
                                     kkt_residual);
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
  // kkt_residual.lq().noalias() += kkt_matrix.Pq().transpose() * s.xi_stack();
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
      ImpulseKKTMatrix& kkt_matrix, ImpulseKKTResidual& kkt_residual) {
  robot.computeImpulseConditionResidual(impulse_status,
                                        impulse_status.contactPoints(),
                                        kkt_residual.P());
  robot.computeImpulseConditionDerivative(impulse_status, kkt_matrix.Pq());
}


inline void ImpulseDynamicsForwardEuler::condenseImpulseDynamics(
    Robot& robot, const ImpulseStatus& impulse_status, 
    ImpulseKKTMatrix& kkt_matrix, ImpulseKKTResidual& kkt_residual) {
  robot.computeMJtJinv(data_.dImDddv, data_.dCdv(), data_.MJtJinv());
  condensing(robot, impulse_status, data_, kkt_matrix, kkt_residual);
}


inline void ImpulseDynamicsForwardEuler::condensing(
    const Robot& robot, const ImpulseStatus& impulse_status, 
    ImpulseDynamicsForwardEulerData& data, ImpulseKKTMatrix& kkt_matrix, 
    ImpulseKKTResidual& kkt_residual) {
  const int dimv = robot.dimv();
  const int dimf = impulse_status.dimp();
  data.MJtJinv_dImDCdqv().leftCols(dimv).noalias() 
      = data.MJtJinv() * data.dImDCdq();
  data.MJtJinv_dImDCdqv().topRightCorner(dimv, dimv).noalias() 
      = data.MJtJinv().topRightCorner(dimv, dimf) * data.dCdv();
  data.MJtJinv_dImDCdqv().bottomRightCorner(dimf, dimv).noalias() 
      = data.MJtJinv().bottomRightCorner(dimf, dimf) * data.dCdv();
  data.MJtJinv_ImDC().noalias() = data.MJtJinv() * data.ImDC();
  data.Qdvfqv().noalias() 
      = (- kkt_matrix.Qdvdvff().diagonal()).asDiagonal() * data.MJtJinv_dImDCdqv();
  data.ldv() = kkt_residual.ldv;
  data.lf() = - kkt_residual.lf();
  data.ldvf().noalias() 
      -= kkt_matrix.Qdvdvff().diagonal().asDiagonal() * data.MJtJinv_ImDC();
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
    const Robot& robot, const ImpulseKKTMatrix& kkt_matrix, 
    const ImpulseKKTResidual& kkt_residual, 
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
    const ImpulseKKTMatrix& kkt_matrix, const ImpulseKKTResidual& kkt_residual, 
    const Eigen::MatrixBase<VectorType>& dgmm, ImpulseSplitDirection& d) {
  assert(dgmm.size() == robot.dimv());
  data.ldvf().noalias() += data.Qdvfqv() * d.dx();
  data.ldv().noalias() += dgmm;
  d.dbetamu().noalias() = - data.MJtJinv() * data.ldvf();
}


inline void ImpulseDynamicsForwardEuler::computeImpulseDynamicsResidual(
    Robot& robot, const ImpulseStatus& impulse_status,
    const ImpulseSplitSolution& s, ImpulseKKTResidual& kkt_residual) {
  setImpulseStatus(impulse_status);
  robot.setImpulseForces(impulse_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, data_.ImD());
  robot.computeImpulseVelocityResidual(impulse_status, data_.C());
  robot.computeImpulseConditionResidual(impulse_status, 
                                        impulse_status.contactPoints(), 
                                        kkt_residual.P());
}


inline double ImpulseDynamicsForwardEuler::l1NormImpulseDynamicsResidual(
    const ImpulseKKTResidual& kkt_residual) const {
  return (data_.ImDC().lpNorm<1>() + kkt_residual.P().lpNorm<1>());
}


inline double ImpulseDynamicsForwardEuler::squaredNormImpulseDynamicsResidual(
    const ImpulseKKTResidual& kkt_residual) const {
  // printf("Error in impulse condition constraints = %lf\n", kkt_residual.P().squaredNorm());
  // return (data_.ImDC().squaredNorm() + kkt_residual.P().squaredNorm());
  return data_.ImDC().squaredNorm();
}


inline void ImpulseDynamicsForwardEuler::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  data_.setImpulseStatus(impulse_status);
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HXX_ 