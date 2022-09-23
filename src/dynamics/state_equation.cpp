#include "robotoc/dynamics/state_equation.hpp"

#include <cassert>


namespace robotoc {

void evalStateEquation(const Robot& robot, const double dt, 
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


void linearizeStateEquation(const Robot& robot, StateEquationData& data,
                            const double dt, const Eigen::VectorXd& q_prev, 
                            const SplitSolution& s, 
                            const SplitSolution& s_next, 
                            SplitKKTMatrix& kkt_matrix, 
                            SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  assert(q_prev.size() == robot.dimq());
  evalStateEquation(robot, dt, s, s_next.q, s_next.v, kkt_residual);
  if (robot.hasFloatingBase()) {
    robot.dSubtractConfiguration_dqf(s.q, s_next.q, kkt_matrix.Fqq());
    data.Fqq_prev.setZero();
    robot.dSubtractConfiguration_dq0(q_prev, s.q, data.Fqq_prev);
    kkt_residual.lq().template head<6>().noalias() 
        += kkt_matrix.Fqq().template topLeftCorner<6, 6>().transpose() 
              * s_next.lmd.template head<6>();
    kkt_residual.lq().template head<6>().noalias() 
        += data.Fqq_prev.template topLeftCorner<6, 6>().transpose() 
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


void correctLinearizeStateEquation(const Robot& robot, StateEquationData& data,
                                   const double dt, const SplitSolution& s, 
                                   const SplitSolution& s_next, 
                                   SplitKKTMatrix& kkt_matrix, 
                                   SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  if (!data.hasFloatingBase()) return;

  data.se3_jac_inverse.compute(data.Fqq_prev, data.Fqq_prev_inv);
  robot.dSubtractConfiguration_dq0(s.q, s_next.q, data.Fqq_prev);
  data.se3_jac_inverse.compute(data.Fqq_prev, data.Fqq_inv);
  data.Fqq_tmp = kkt_matrix.Fqq().template topLeftCorner<6, 6>();
  kkt_matrix.Fqq().template topLeftCorner<6, 6>().noalias() = - data.Fqq_inv * data.Fqq_tmp;
  kkt_matrix.Fqv().template topLeftCorner<6, 6>() = - dt * data.Fqq_inv;
  data.Fq_tmp  = kkt_residual.Fq().template head<6>();
  kkt_residual.Fq().template head<6>().noalias() = - data.Fqq_inv * data.Fq_tmp;
  data.Fq_tmp = kkt_matrix.fq().template head<6>();
  kkt_matrix.fq().template head<6>().noalias() = - data.Fqq_inv * data.Fq_tmp;
}


void correctCostateDirection(StateEquationData& data, SplitDirection& d)  {
  if (!data.hasFloatingBase()) return;

  data.Fq_tmp.noalias() = data.Fqq_prev_inv.transpose() * d.dlmdgmm.template head<6>();
  d.dlmdgmm.template head<6>() = - data.Fq_tmp;
}


void computeInitialStateDirection(const Robot& robot, 
                                  const StateEquationData& data, 
                                  const Eigen::VectorXd& q0, 
                                  const Eigen::VectorXd& v0, 
                                  const SplitSolution& s0, SplitDirection& d0) {
  robot.subtractConfiguration(q0, s0.q, d0.dq());
  if (robot.hasFloatingBase()) {
    d0.dq().template head<6>().noalias() 
        = - data.Fqq_prev_inv * d0.dq().template head<6>();
  }
  d0.dv() = v0 - s0.v;
}


StateEquation::StateEquation(const Robot& robot)
  : data_(robot) {
}


StateEquation::StateEquation()
  : data_() {
}



void StateEquation::evalStateEquation(const Robot& robot, const double dt, 
                                      const SplitSolution& s, 
                                      const Eigen::VectorXd& q_next, 
                                      const Eigen::VectorXd& v_next, 
                                      SplitKKTResidual& kkt_residual) const {
  ::robotoc::evalStateEquation(robot, dt, s, q_next, v_next, kkt_residual);
}


void StateEquation::linearizeStateEquation(const Robot& robot, const double dt, 
                                           const Eigen::VectorXd& q_prev, 
                                           const SplitSolution& s, 
                                           const SplitSolution& s_next, 
                                           SplitKKTMatrix& kkt_matrix, 
                                           SplitKKTResidual& kkt_residual) {
  ::robotoc::linearizeStateEquation(robot, data_, dt, q_prev, s, s_next, kkt_matrix, kkt_residual);
}


void StateEquation::correctLinearizedStateEquation(
    const Robot& robot, const double dt, const SplitSolution& s, 
    const SplitSolution& s_next, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) {
  ::robotoc::correctLinearizeStateEquation(robot, data_, dt, s, s_next, kkt_matrix, kkt_residual);
}


void StateEquation::correctCostateDirection(SplitDirection& d)  {
  ::robotoc::correctCostateDirection(data_, d);
}


void StateEquation::computeInitialStateDirection(const Robot& robot, 
                                                 const Eigen::VectorXd& q0, 
                                                 const Eigen::VectorXd& v0, 
                                                 const SplitSolution& s0, 
                                                 SplitDirection& d0) const {
  ::robotoc::computeInitialStateDirection(robot, data_, q0, v0, s0, d0);
}

} // namespace robotoc 