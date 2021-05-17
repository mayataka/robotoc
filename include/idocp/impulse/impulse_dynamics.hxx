#ifndef IDOCP_IMPULSE_DYNAMICS_HXX_
#define IDOCP_IMPULSE_DYNAMICS_HXX_

#include "idocp/impulse/impulse_dynamics.hpp"

#include <cassert>

namespace idocp {

inline ImpulseDynamics::ImpulseDynamics(
    const Robot& robot) 
  : data_(robot) {
}


inline ImpulseDynamics::ImpulseDynamics() 
  : data_() {
}


inline ImpulseDynamics::~ImpulseDynamics() {
}


inline void ImpulseDynamics::linearizeImpulseDynamics(
    Robot& robot, const ImpulseStatus& impulse_status,  
    const ImpulseSplitSolution& s, ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual) {
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
}


inline void ImpulseDynamics::linearizeInverseImpulseDynamics(
    Robot& robot, const ImpulseStatus& impulse_status, 
    const ImpulseSplitSolution& s, ImpulseDynamicsData& data) {
  robot.setImpulseForces(impulse_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, data.ImD());
  robot.RNEAImpulseDerivatives(s.q, s.dv, data.dImDdq(), data.dImDddv);
}


inline void ImpulseDynamics::linearizeImpulseVelocityConstraint(
    Robot& robot, const ImpulseStatus& impulse_status, 
    ImpulseDynamicsData& data) {
  robot.computeImpulseVelocityResidual(impulse_status, data.C());
  robot.computeImpulseVelocityDerivatives(impulse_status, data.dCdq(), 
                                          data.dCdv());
}


inline void ImpulseDynamics::condenseImpulseDynamics(
    Robot& robot, const ImpulseStatus& impulse_status, 
    ImpulseSplitKKTMatrix& kkt_matrix, ImpulseSplitKKTResidual& kkt_residual) {
  robot.computeMJtJinv(data_.dImDddv, data_.dCdv(), data_.MJtJinv());
  condensing(robot, impulse_status, data_, kkt_matrix, kkt_residual);
}


inline void ImpulseDynamics::condensing(
    const Robot& robot, const ImpulseStatus& impulse_status, 
    ImpulseDynamicsData& data, ImpulseSplitKKTMatrix& kkt_matrix, 
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
  data.Qdvfqv().bottomLeftCorner(dimf, dimv).noalias() 
      -= kkt_matrix.Qqf().transpose();

  data.ldv() = kkt_residual.ldv;
  data.lf()  = - kkt_residual.lf();
  data.ldv().noalias() 
      -= kkt_matrix.Qdvdv().diagonal().asDiagonal() 
          * data.MJtJinv_ImDC().head(dimv);
  data.lf().noalias() -= kkt_matrix.Qff() * data.MJtJinv_ImDC().tail(dimf);

  kkt_matrix.Qxx().noalias() 
      -= data.MJtJinv_dImDCdqv().transpose() * data.Qdvfqv();
  kkt_matrix.Qxx().topRows(dimv).noalias() 
      += kkt_matrix.Qqf() * data.MJtJinv_dImDCdqv().bottomRows(dimf);

  kkt_residual.lx().noalias() 
      -= data.MJtJinv_dImDCdqv().transpose() * data.ldvf();
  kkt_residual.lq().noalias()
      += kkt_matrix.Qqf() * data.MJtJinv_ImDC().tail(dimf);

  kkt_matrix.Fvq() = - data.MJtJinv_dImDCdqv().topLeftCorner(dimv, dimv);
  kkt_matrix.Fvv() = Eigen::MatrixXd::Identity(dimv, dimv) 
                    - data.MJtJinv_dImDCdqv().topRightCorner(dimv, dimv);
  kkt_residual.Fv().noalias() -= data.MJtJinv_ImDC().head(dimv);
}


inline void ImpulseDynamics::computeCondensedPrimalDirection(
    const Robot& robot, ImpulseSplitDirection& d) {
  expansionPrimal(robot, data_, d);
}


template <typename VectorType>
inline void ImpulseDynamics::computeCondensedDualDirection(
    const Robot& robot, const ImpulseSplitKKTMatrix& kkt_matrix, 
    const ImpulseSplitKKTResidual& kkt_residual, 
    const Eigen::MatrixBase<VectorType>& dgmm, ImpulseSplitDirection& d) {
  assert(dgmm.size() == robot.dimv());
  expansionDual(robot, data_, kkt_matrix, kkt_residual, dgmm, d);
}


inline void ImpulseDynamics::expansionPrimal(
    const Robot& robot, const ImpulseDynamicsData& data, 
    ImpulseSplitDirection& d) {
  d.ddvf().noalias()  = - data.MJtJinv_dImDCdqv() * d.dx();
  d.ddvf().noalias() -= data.MJtJinv_ImDC();
  d.df().array() *= -1;
}


template <typename VectorType>
inline void ImpulseDynamics::expansionDual(
    const Robot& robot, ImpulseDynamicsData& data, 
    const ImpulseSplitKKTMatrix& kkt_matrix, 
    const ImpulseSplitKKTResidual& kkt_residual, 
    const Eigen::MatrixBase<VectorType>& dgmm, ImpulseSplitDirection& d) {
  assert(dgmm.size() == robot.dimv());
  data.ldvf().noalias() += data.Qdvfqv() * d.dx();
  data.ldv().noalias()  += dgmm;
  d.dbetamu().noalias() = - data.MJtJinv() * data.ldvf();
}


inline void ImpulseDynamics::computeImpulseDynamicsResidual(
    Robot& robot, const ImpulseStatus& impulse_status,
    const ImpulseSplitSolution& s, ImpulseSplitKKTResidual& kkt_residual) {
  setImpulseStatus(impulse_status);
  robot.setImpulseForces(impulse_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, data_.ImD());
  robot.computeImpulseVelocityResidual(impulse_status, data_.C());
}


inline double ImpulseDynamics::l1NormImpulseDynamicsResidual(
    const ImpulseSplitKKTResidual& kkt_residual) const {
  return data_.ImDC().lpNorm<1>();
}


inline double ImpulseDynamics::squaredNormImpulseDynamicsResidual(
    const ImpulseSplitKKTResidual& kkt_residual) const {
  return data_.ImDC().squaredNorm();
}


inline void ImpulseDynamics::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  data_.setImpulseStatus(impulse_status);
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_DYNAMICS_HXX_ 