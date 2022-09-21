#include "robotoc/ocp/state_equation.hpp"

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


StateEquation::~StateEquation() {
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
  linearizeStateEquation_impl(robot, dt, q_prev, s, s_next, kkt_matrix, kkt_residual);
}


void StateEquation::correctLinearizedStateEquation(
    const Robot& robot, const double dt, const SplitSolution& s, 
    const SplitSolution& s_next, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) {
  correctLinearizedStateEquation_impl(robot, dt, s, s_next, 
                                      kkt_matrix, kkt_residual);
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