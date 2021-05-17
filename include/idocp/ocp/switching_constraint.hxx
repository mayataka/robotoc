#ifndef IDOCP_SWITCHING_CONSTRAINT_HXX_ 
#define IDOCP_SWITCHING_CONSTRAINT_HXX_

#include "idocp/ocp/switching_constraint.hpp"

#include <cassert>

namespace idocp {

inline SwitchingConstraint::SwitchingConstraint(
    const Robot& robot) 
  : q_(Eigen::VectorXd::Zero(robot.dimq())), 
    dq_(Eigen::VectorXd::Zero(robot.dimv())) {
}


inline SwitchingConstraint::SwitchingConstraint() 
  : q_(), 
    dq_() {
}


inline SwitchingConstraint::~SwitchingConstraint() {
}


inline void SwitchingConstraint::linearizeSwitchingConstraint(
    Robot& robot, const ImpulseStatus& impulse_status, const double dt1, 
    const double dt2, const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual, SplitStateConstraintJacobian& jac) {
  assert(dt1 > 0);
  assert(dt2 > 0);
  jac.setImpulseStatus(impulse_status);
  computeSwitchingConstraintResidual(robot, impulse_status, dt1, dt2, s, 
                                     kkt_residual);
  robot.computeContactPositionDerivative(impulse_status, kkt_matrix.Pq());
  if (robot.hasFloatingBase()) {
    robot.dIntegratedConfiguration(s.q, dq_, jac.dintegrate_dq);
    robot.dIntegratedVelocity(s.q, dq_, jac.dintegrate_dv);
    jac.Phiq().noalias() = kkt_matrix.Pq() * jac.dintegrate_dq;
    jac.Phiv().noalias() = (dt1+dt2) * kkt_matrix.Pq() * jac.dintegrate_dv;
    jac.Phia().noalias() = (dt1*dt2) * kkt_matrix.Pq() * jac.dintegrate_dv;
  }
  else {
    jac.Phiq() = kkt_matrix.Pq();
    jac.Phiv() = (dt1+dt2) * kkt_matrix.Pq();
    jac.Phia() = (dt1*dt2) * kkt_matrix.Pq();
  }
  kkt_residual.lq().noalias() += jac.Phiq().transpose() * s.xi_stack();
  kkt_residual.lv().noalias() += jac.Phiv().transpose() * s.xi_stack();
  kkt_residual.la.noalias()   += jac.Phia().transpose() * s.xi_stack();
}


inline void SwitchingConstraint::computeSwitchingConstraintResidual(
    Robot& robot, const ImpulseStatus& impulse_status, const double dt1, 
    const double dt2, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) {
  assert(dt1 > 0);
  assert(dt2 > 0);
  dq_ = (dt1+dt2) * s.v + (dt1*dt2) * s.a;
  robot.integrateConfiguration(s.q, dq_, 1.0, q_);
  robot.updateKinematics(q_);
  robot.computeContactPositionResidual(impulse_status, 
                                       impulse_status.contactPoints(), 
                                       kkt_residual.P());
}


inline double SwitchingConstraint::l1NormSwitchingConstraintResidual(
    const SplitKKTResidual& kkt_residual) {
  return kkt_residual.P().template lpNorm<1>();
}


inline double SwitchingConstraint::squaredNormSwitchingConstraintResidual(
    const SplitKKTResidual& kkt_residual) {
  return kkt_residual.P().squaredNorm();
}

} // namespace idocp

#endif // IDOCP_SWITCHING_CONSTRAINT_HXX_