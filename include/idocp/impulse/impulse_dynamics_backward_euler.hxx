#ifndef IDOCP_IMPULSE_DYNAMICS_BACKWARD_EULER_HXX
#define IDOCP_IMPULSE_DYNAMICS_BACKWARD_EULER_HXX

#include "idocp/impulse/impulse_dynamics_backward_euler.hpp"

#include <assert.h>

namespace idocp {

inline ImpulseDynamicsBackwardEuler::ImpulseDynamicsBackwardEuler(
    const Robot& robot) 
  : data_(robot),
    dimv_(robot.dimv()),
    dimf_(0) {
}


inline ImpulseDynamicsBackwardEuler::ImpulseDynamicsBackwardEuler() 
  : data_(), 
    dimv_(0), 
    dimf_(0) {
}


inline ImpulseDynamicsBackwardEuler::~ImpulseDynamicsBackwardEuler() {
}


inline void ImpulseDynamicsBackwardEuler::linearizeImpulseDynamics(
    Robot& robot, const ImpulseStatus& impulse_status,  
    const ImpulseSplitSolution& s, ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual) {
  setImpulseStatus(impulse_status);
  linearizeInverseImpulseDynamics(robot, impulse_status, s, data_);
  linearizeImpulseVelocityConstraint(robot, impulse_status, kkt_matrix, 
                                     kkt_residual);
  // augment inverse impulse dynamics constraint
  kkt_residual.lq().noalias() += data_.dImDdq.transpose() * s.beta;
  kkt_residual.ldv.noalias()  += data_.dImDddv.transpose() * s.beta;
  // We use an equivalence dmIDdf().transpose() = - Vv() to avoid
  // redundant calculation of dImDdf().
  kkt_residual.lf().noalias() -= kkt_matrix.Vv() * s.beta;
  // augment impulse velocity constraint
  kkt_residual.lq().noalias() += kkt_matrix.Vq().transpose() * s.mu_stack();
  kkt_residual.lv().noalias() += kkt_matrix.Vv().transpose() * s.mu_stack();
}


inline void ImpulseDynamicsBackwardEuler::linearizeInverseImpulseDynamics(
    Robot& robot, const ImpulseStatus& impulse_status, 
    const ImpulseSplitSolution& s, ImpulseDynamicsBackwardEulerData& data) {
  robot.setImpulseForces(impulse_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, data.ImD);
  robot.RNEAImpulseDerivatives(s.q, s.dv, data.dImDdq, data.dImDddv);
}


inline void ImpulseDynamicsBackwardEuler::linearizeImpulseVelocityConstraint(
    Robot& robot, const ImpulseStatus& impulse_status, 
    ImpulseSplitKKTMatrix& kkt_matrix, ImpulseSplitKKTResidual& kkt_residual) {
  robot.computeImpulseVelocityResidual(impulse_status, kkt_residual.V());
  robot.computeImpulseVelocityDerivatives(impulse_status, kkt_matrix.Vq(), 
                                          kkt_matrix.Vv());
}


inline void ImpulseDynamicsBackwardEuler::condenseImpulseDynamics(
    Robot& robot, const ImpulseStatus& impulse_status, 
    ImpulseSplitKKTMatrix& kkt_matrix, ImpulseSplitKKTResidual& kkt_residual) {
  robot.computeMinv(data_.dImDddv, data_.Minv);
  condensing(robot, data_, kkt_matrix, kkt_residual);
}


inline void ImpulseDynamicsBackwardEuler::condensing(
    const Robot& robot, ImpulseDynamicsBackwardEulerData& data, 
    ImpulseSplitKKTMatrix& kkt_matrix, ImpulseSplitKKTResidual& kkt_residual) {
  const int dimv = robot.dimv();
  kkt_matrix.Fvq().noalias() = data.Minv * data.dImDdq;
  kkt_matrix.Fvf().noalias() = - data.Minv * kkt_matrix.Vv().transpose(); // this is Minv_dImDdf
  data.Minv_ImD.noalias() = data.Minv * data.ImD;
  data.Qdvq.noalias() 
    = (- kkt_matrix.Qdvdv().diagonal()).asDiagonal() * kkt_matrix.Fvq();
  data.Qdvf().noalias() 
    = (- kkt_matrix.Qdvdv().diagonal()).asDiagonal() * kkt_matrix.Fvf();
  data.ldv.noalias() 
    = (- kkt_matrix.Qdvdv().diagonal()).asDiagonal() * data.Minv_ImD;
  kkt_matrix.Qqq().noalias() -= kkt_matrix.Fvq().transpose() * data.Qdvq;
  kkt_matrix.Qfq().transpose().noalias() 
        -= kkt_matrix.Fvq().transpose() * data.Qdvf();
  kkt_matrix.Qff().noalias() -= kkt_matrix.Fvf().transpose() * data.Qdvf();
  kkt_residual.lq().noalias() -= kkt_matrix.Fvq().transpose() * data.ldv;
  kkt_residual.lf().noalias() -= kkt_matrix.Fvf().transpose() * data.ldv;
  kkt_matrix.Fvv() = - Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_residual.Fv().noalias() -= data.Minv_ImD;
}


inline void ImpulseDynamicsBackwardEuler::computeCondensedPrimalDirection(
    const Robot& robot, const ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitDirection& d) const {
  d.ddv() = - data_.Minv_ImD;
  d.ddv().noalias() -= kkt_matrix.Fvq() * d.dq();
  d.ddv().noalias() -= kkt_matrix.Fvf() * d.df();
}


inline void ImpulseDynamicsBackwardEuler::computeCondensedDualDirection(
    const Robot& robot, ImpulseSplitDirection& d) {
  data_.ldv.noalias() += data_.Qdvq * d.dq();
  data_.ldv.noalias() += data_.Qdvf() * d.df();
  data_.ldv.noalias() += d.dgmm();
  d.dbeta().noalias() = - data_.Minv * data_.ldv;
}


inline void ImpulseDynamicsBackwardEuler::computeImpulseDynamicsResidual(
    Robot& robot, const ImpulseStatus& impulse_status,
    const ImpulseSplitSolution& s, ImpulseSplitKKTResidual& kkt_residual) {
  setImpulseStatus(impulse_status);
  robot.setImpulseForces(impulse_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, data_.ImD);
  robot.computeImpulseVelocityResidual(impulse_status, kkt_residual.V());
}


inline double ImpulseDynamicsBackwardEuler::l1NormImpulseDynamicsResidual(
    const ImpulseSplitKKTResidual& kkt_residual) const {
  return (data_.ImD.lpNorm<1>() + kkt_residual.V().lpNorm<1>());
}


inline double ImpulseDynamicsBackwardEuler::squaredNormImpulseDynamicsResidual(
    const ImpulseSplitKKTResidual& kkt_residual) const {
  return (data_.ImD.squaredNorm() + kkt_residual.V().squaredNorm());
}


inline void ImpulseDynamicsBackwardEuler::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  data_.setImpulseStatus(impulse_status);
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_DYNAMICS_BACKWARD_EULER_HXX 