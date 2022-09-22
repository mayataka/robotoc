#include "robotoc/dynamics/terminal_state_equation.hpp"

#include <cassert>

namespace robotoc {

TerminalStateEquation::TerminalStateEquation(const Robot& robot)
  : data_(robot),
    has_floating_base_(robot.hasFloatingBase()) {
}


TerminalStateEquation::TerminalStateEquation()
  : data_(),
    has_floating_base_(false) {
}


void TerminalStateEquation::linearizeStateEquation(
    const Robot& robot, const Eigen::VectorXd& q_prev, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  if (robot.hasFloatingBase()) {
    data_.Fqq_prev.setZero();
    robot.dSubtractConfiguration_dq0(q_prev, s.q, data_.Fqq_prev);
    kkt_residual.lq().template head<6>().noalias() 
        += data_.Fqq_prev.template topLeftCorner<6, 6>().transpose() 
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
  if (!has_floating_base_) return;

  data_.se3_jac_inverse.compute(data_.Fqq_prev, data_.Fqq_prev_inv);
}


void TerminalStateEquation::correctCostateDirection(SplitDirection& d) {
  if (!has_floating_base_) return;

  data_.Fq_tmp = data_.Fqq_prev_inv.transpose() * d.dlmdgmm.template head<6>();
  d.dlmdgmm.template head<6>() = - data_.Fq_tmp;
}

} // namespace robotoc 