#include "robotoc/dynamics/unconstr_dynamics.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

UnconstrDynamics::UnconstrDynamics(const Robot& robot) 
  : lu_condensed_(Eigen::VectorXd::Zero(robot.dimv())),
    ID_(Eigen::VectorXd::Zero(robot.dimv())),
    Quu_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dID_dq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dID_dv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dID_da_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Quu_dID_dq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Quu_dID_dv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Quu_dID_da_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dimv_(robot.dimv()) {
  if (robot.hasFloatingBase()) {
    throw std::invalid_argument(
        "[UnconstrDynamics] robot has a floating base: robot should have no constraints!");
  }
  if (robot.maxNumContacts() > 0) {
    throw std::invalid_argument(
        "[UnconstrDynamics] robot have contact frames: robot should have no constraints!");
  }
}


UnconstrDynamics::UnconstrDynamics() 
  : lu_condensed_(),
    ID_(),
    Quu_(),
    dID_dq_(),
    dID_dv_(),
    dID_da_(),
    Quu_dID_dq_(),
    Quu_dID_dv_(),
    Quu_dID_da_(),
    dimv_(0) {
}


void UnconstrDynamics::evalUnconstrDynamics(Robot& robot, 
                                            const SplitSolution& s) {
  robot.RNEA(s.q, s.v, s.a, ID_);
  ID_.noalias() -= s.u;
}


void UnconstrDynamics::linearizeUnconstrDynamics(Robot& robot, const double dt, 
                                                 const SplitSolution& s, 
                                                 SplitKKTResidual& kkt_residual) { 
  assert(dt > 0);
  evalUnconstrDynamics(robot, s);
  // augment inverse dynamics constraint
  robot.RNEADerivatives(s.q, s.v, s.a, dID_dq_, dID_dv_, dID_da_);
  kkt_residual.lq().noalias() += dt * dID_dq_.transpose() * s.beta;
  kkt_residual.lv().noalias() += dt * dID_dv_.transpose() * s.beta;
  kkt_residual.la.noalias()   += dt * dID_da_.transpose() * s.beta;
  kkt_residual.lu.noalias()   -= dt * s.beta; 
}


void UnconstrDynamics::condenseUnconstrDynamics(SplitKKTMatrix& kkt_matrix, 
                                                SplitKKTResidual& kkt_residual) {
  // condense KKT residual
  lu_condensed_ = kkt_residual.lu;
  lu_condensed_.noalias() += kkt_matrix.Quu.diagonal().asDiagonal() * ID_;
  kkt_residual.lq().noalias() += dID_dq_.transpose() * lu_condensed_;
  kkt_residual.lv().noalias() += dID_dv_.transpose() * lu_condensed_;
  kkt_residual.la.noalias()   += dID_da_.transpose() * lu_condensed_;
  // condense KKT Hessian
  Quu_dID_dq_.noalias() = kkt_matrix.Quu.diagonal().asDiagonal() * dID_dq_;
  Quu_dID_dv_.noalias() = kkt_matrix.Quu.diagonal().asDiagonal() * dID_dv_;
  Quu_dID_da_.noalias() = kkt_matrix.Quu.diagonal().asDiagonal() * dID_da_;
  kkt_matrix.Qqq().noalias() += dID_dq_.transpose() * Quu_dID_dq_;
  kkt_matrix.Qqv().noalias() += dID_dq_.transpose() * Quu_dID_dv_;
  kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
  kkt_matrix.Qvv().noalias() += dID_dv_.transpose() * Quu_dID_dv_;
  kkt_matrix.Qaa.noalias()   += dID_da_.transpose() * Quu_dID_da_;
  //  this is actually Qqa and Qqv
  kkt_matrix.Qqu().noalias() = Quu_dID_dq_.transpose() * dID_da_;
  kkt_matrix.Qvu().noalias() = Quu_dID_dv_.transpose() * dID_da_;
}


void UnconstrDynamics::expandPrimal(SplitDirection& d) const {
  d.du = ID_;
  d.du.noalias() += dID_dq_ * d.dq();
  d.du.noalias() += dID_dv_ * d.dv();
  d.du.noalias() += dID_da_ * d.da();
}


void UnconstrDynamics::expandDual(const double dt, 
                                  const SplitKKTMatrix& kkt_matrix, 
                                  const SplitKKTResidual& kkt_residual, 
                                  SplitDirection& d) {
  assert(dt > 0);
  d.dbeta().noalias() = (kkt_residual.lu  + kkt_matrix.Quu * d.du) / dt;
}

} // namespace robotoc 