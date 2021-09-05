#ifndef IDOCP_TERMINAL_STATE_EQUATION_HXX_
#define IDOCP_TERMINAL_STATE_EQUATION_HXX_

#include "idocp/ocp/terminal_state_equation.hpp"

#include <cassert>

namespace idocp {

inline TerminalStateEquation::TerminalStateEquation(const Robot& robot)
  : Fqq_prev_inv_(),
    Fq_tmp_(),
    lie_der_inverter_(),
    has_floating_base_(robot.hasFloatingBase()) {
  if (robot.hasFloatingBase()) {
    Fqq_prev_inv_.resize(6, 6);
    Fqq_prev_inv_.setZero();
    Fq_tmp_.resize(6);
    Fq_tmp_.setZero();
  }
}


inline TerminalStateEquation::TerminalStateEquation()
  : Fqq_prev_inv_(),
    Fq_tmp_(),
    lie_der_inverter_(),
    has_floating_base_(false) {
}


inline TerminalStateEquation::~TerminalStateEquation() {
}


template <typename ConfigVectorType>
inline void TerminalStateEquation::linearizeStateEquation(
    const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
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


template <typename ConfigVectorType>
inline void TerminalStateEquation::linearizeStateEquationAlongLieGroup(
    const Robot& robot, const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) {
  linearizeStateEquation(robot, q_prev, s, kkt_matrix, kkt_residual);
  if (has_floating_base_) {
    lie_der_inverter_.computeLieDerivativeInverse(kkt_matrix.Fqq_prev, 
                                                  Fqq_prev_inv_);
  }
}


inline void TerminalStateEquation::correctCostateDirection(SplitDirection& d) {
  if (has_floating_base_) {
    Fq_tmp_ = Fqq_prev_inv_.transpose() * d.dlmdgmm.template head<6>();
    d.dlmdgmm.template head<6>() = - Fq_tmp_;
  }
}

} // namespace idocp 

#endif // IDOCP_TERMINAL_STATE_EQUATION_HXX_ 