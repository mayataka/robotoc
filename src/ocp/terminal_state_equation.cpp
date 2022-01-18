#include "robotoc/ocp/terminal_state_equation.hpp"

#include <cassert>

namespace robotoc {

TerminalStateEquation::TerminalStateEquation(const Robot& robot)
  : Fqq_prev_inv_(),
    Fq_tmp_(),
    se3_jac_inverse_(),
    has_floating_base_(robot.hasFloatingBase()) {
  if (robot.hasFloatingBase()) {
    Fqq_prev_inv_.resize(6, 6);
    Fqq_prev_inv_.setZero();
    Fq_tmp_.resize(6);
    Fq_tmp_.setZero();
  }
}


TerminalStateEquation::TerminalStateEquation()
  : Fqq_prev_inv_(),
    Fq_tmp_(),
    se3_jac_inverse_(),
    has_floating_base_(false) {
}


TerminalStateEquation::~TerminalStateEquation() {
}


void TerminalStateEquation::linearizeStateEquation(
    const Robot& robot, const Eigen::VectorXd& q_prev, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  if (robot.hasFloatingBase()) {
    robot.dSubtractConfiguration_dq0(q_prev, s.q, kkt_matrix.Fqq_prev);
    kkt_residual.lq().template head<6>().noalias() 
        += kkt_matrix.Fqq_prev.template topLeftCorner<6, 6>().transpose() 
              * s.lmd.template head<6>();
    kkt_residual.lq().tail(robot.dimv()-6).noalias() 
        -= s.lmd.tail(robot.dimv()-6);
  }
  else {
    kkt_residual.lq().noalias() -= s.lmd;
  }
  kkt_residual.lv().noalias() -= s.gmm;
}


void TerminalStateEquation::correctLinearizedStateEquation(
    SplitKKTMatrix& kkt_matrix) {
  if (has_floating_base_) {
    se3_jac_inverse_.compute(kkt_matrix.Fqq_prev, Fqq_prev_inv_);
  }
}


void TerminalStateEquation::correctCostateDirection(SplitDirection& d) {
  if (has_floating_base_) {
    Fq_tmp_ = Fqq_prev_inv_.transpose() * d.dlmdgmm.template head<6>();
    d.dlmdgmm.template head<6>() = - Fq_tmp_;
  }
}

} // namespace robotoc 