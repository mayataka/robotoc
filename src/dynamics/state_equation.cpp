#include "robotoc/dynamics/state_equation.hpp"

#include <cassert>


namespace robotoc {

StateEquation::StateEquation(const Robot& robot)
  : Fqq_inv_(),
    Fqq_prev_inv_(),
    Fqq_tmp_(),
    Fq_tmp_(),
    se3_jac_inverse_(),
    has_floating_base_(robot.hasFloatingBase()) {
  if (robot.hasFloatingBase()) {
    Fqq_inv_.resize(6, 6);
    Fqq_inv_.setZero();
    Fqq_prev_inv_.resize(6, 6);
    Fqq_prev_inv_.setZero();
    Fqq_tmp_.resize(6, 6);
    Fqq_tmp_.setZero();
    Fq_tmp_.resize(6);
    Fq_tmp_.setZero();
  }
}


StateEquation::StateEquation()
  : Fqq_inv_(),
    Fqq_prev_inv_(),
    Fqq_tmp_(),
    Fq_tmp_(),
    se3_jac_inverse_(),
    has_floating_base_(false) {
}


void StateEquation::evalStateEquation(const Robot& robot, const double dt, 
                                      const SplitSolution& s, 
                                      const Eigen::VectorXd& q_next, 
                                      const Eigen::VectorXd& v_next, 
                                      SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  assert(q_next.size() == robot.dimq());
  assert(v_next.size() == robot.dimv());
  robot.subtractConfiguration(s.q, q_next, kkt_residual.Fq());
  kkt_residual.Fq().noalias() += dt * s.v;
  kkt_residual.Fv() = s.v + dt * s.a - v_next;
}


void StateEquation::linearizeStateEquation(const Robot& robot, const double dt, 
                                           const Eigen::VectorXd& q_prev, 
                                           const SplitSolution& s, 
                                           const SplitSolution& s_next, 
                                           SplitKKTMatrix& kkt_matrix, 
                                           SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  assert(q_prev.size() == robot.dimq());
  evalStateEquation(robot, dt, s, s_next.q, s_next.v, kkt_residual);
  if (robot.hasFloatingBase()) {
    robot.dSubtractConfiguration_dqf(s.q, s_next.q, kkt_matrix.Fqq());
    robot.dSubtractConfiguration_dq0(q_prev, s.q, kkt_matrix.Fqq_prev);
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
  kkt_matrix.Fqv().diagonal().fill(dt);
  kkt_residual.lv().noalias() += dt * s_next.lmd + s_next.gmm - s.gmm;
  kkt_residual.la.noalias() += dt * s_next.gmm;
  // STO sensitivities
  kkt_residual.h += s_next.lmd.dot(s.v);
  kkt_residual.h += s_next.gmm.dot(s.a);
  kkt_matrix.hv().noalias() += s_next.lmd;
  kkt_matrix.ha.noalias()   += s_next.gmm;
  kkt_matrix.fq() = s.v;
  kkt_matrix.fv() = s.a;
}


void StateEquation::correctLinearizedStateEquation(
    const Robot& robot, const double dt, const SplitSolution& s, 
    const SplitSolution& s_next, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  if (!has_floating_base_) return;

  se3_jac_inverse_.compute(kkt_matrix.Fqq_prev, Fqq_prev_inv_);
  robot.dSubtractConfiguration_dq0(s.q, s_next.q, kkt_matrix.Fqq_prev);
  se3_jac_inverse_.compute(kkt_matrix.Fqq_prev, Fqq_inv_);
  Fqq_tmp_ = kkt_matrix.Fqq().template topLeftCorner<6, 6>();
  kkt_matrix.Fqq().template topLeftCorner<6, 6>().noalias() = - Fqq_inv_ * Fqq_tmp_;
  kkt_matrix.Fqv().template topLeftCorner<6, 6>() = - dt * Fqq_inv_;
  Fq_tmp_  = kkt_residual.Fq().template head<6>();
  kkt_residual.Fq().template head<6>().noalias() = - Fqq_inv_ * Fq_tmp_;
  Fq_tmp_ = kkt_matrix.fq().template head<6>();
  kkt_matrix.fq().template head<6>().noalias() = - Fqq_inv_ * Fq_tmp_;
}


void StateEquation::correctCostateDirection(SplitDirection& d)  {
  if (!has_floating_base_) return;

  Fq_tmp_.noalias() = Fqq_prev_inv_.transpose() * d.dlmdgmm.template head<6>();
  d.dlmdgmm.template head<6>() = - Fq_tmp_;
}


void StateEquation::computeInitialStateDirection(const Robot& robot, 
                                                 const Eigen::VectorXd& q0, 
                                                 const Eigen::VectorXd& v0, 
                                                 const SplitSolution& s0, 
                                                 SplitDirection& d0) const {
  robot.subtractConfiguration(q0, s0.q, d0.dq());
  if (robot.hasFloatingBase()) {
    d0.dq().template head<6>().noalias() 
        = - Fqq_prev_inv_ * d0.dq().template head<6>();
  }
  d0.dv() = v0 - s0.v;
}

} // namespace robotoc 