#include "robotoc/dynamics/terminal_state_equation.hpp"

#include <cassert>

namespace robotoc {

TerminalStateEquation::TerminalStateEquation(const Robot& robot)
  : data_(robot) {
}


TerminalStateEquation::TerminalStateEquation()
  : data_() {
}


void TerminalStateEquation::linearize(const Robot& robot, 
                                      StateEquationData& data, 
                                      const Eigen::VectorXd& q_prev, 
                                      const SplitSolution& s, 
                                      SplitKKTMatrix& kkt_matrix, 
                                      SplitKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  if (robot.hasFloatingBase()) {
    data.Fqq_prev.setZero();
    robot.dSubtractConfiguration_dq0(q_prev, s.q, data.Fqq_prev);
    kkt_residual.lq().template head<6>().noalias() 
        += data.Fqq_prev.template topLeftCorner<6, 6>().transpose() 
              * s.lmd.template head<6>();
    kkt_residual.lq().tail(robot.dimv()-6).noalias() 
        -= s.lmd.tail(robot.dimv()-6);
  }
  else {
    kkt_residual.lq().noalias() -= s.lmd;
  }
  kkt_residual.lv().noalias() -= s.gmm;
}


void TerminalStateEquation::correctLinearize(StateEquationData& data, 
                                             SplitKKTMatrix& kkt_matrix) {
  if (!data.hasFloatingBase()) return;

  data.se3_jac_inverse.compute(data.Fqq_prev, data.Fqq_prev_inv);
}


void TerminalStateEquation::correctCostateDirection(StateEquationData& data, 
                                                    SplitDirection& d) {
  if (!data.hasFloatingBase()) return;

  data.Fq_tmp = data.Fqq_prev_inv.transpose() * d.dlmdgmm.template head<6>();
  d.dlmdgmm.template head<6>() = - data.Fq_tmp;
}


void TerminalStateEquation::linearizeStateEquation(
    const Robot& robot, const Eigen::VectorXd& q_prev, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) {
  linearize(robot, data_, q_prev, s, kkt_matrix, kkt_residual);
}


void TerminalStateEquation::correctLinearizedStateEquation(
    SplitKKTMatrix& kkt_matrix) {
  correctLinearize(data_, kkt_matrix);
}


void TerminalStateEquation::correctCostateDirection(SplitDirection& d) {
  correctCostateDirection(data_, d);
}

} // namespace robotoc 