#include "robotoc/dynamics/terminal_state_equation.hpp"
#include "robotoc/dynamics/state_equation.hpp"

#include <cassert>

namespace robotoc {

void linearizeTerminalStateEquation(const Robot& robot, 
                                    const Eigen::VectorXd& q_prev, 
                                    const SplitSolution& s, 
                                    StateEquationData& data, 
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


void correctLinearizeTerminalStateEquation(StateEquationData& data, 
                                           SplitKKTMatrix& kkt_matrix) {
  if (!data.hasFloatingBase()) return;

  data.se3_jac_inverse.compute(data.Fqq_prev, data.Fqq_prev_inv);
}

} // namespace robotoc 