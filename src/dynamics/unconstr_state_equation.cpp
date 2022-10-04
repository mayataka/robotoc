#include "robotoc/dynamics/unconstr_state_equation.hpp"

#include <cassert>


namespace robotoc {

void linearizeUnconstrForwardEuler(const double dt, const SplitSolution& s, 
                                   const SplitSolution& s_next, 
                                   SplitKKTMatrix& kkt_matrix, 
                                   SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  evalUnconstrForwardEuler(dt, s, s_next, kkt_residual);
  kkt_residual.lq().noalias() += s_next.lmd - s.lmd;
  kkt_residual.lv().noalias() += dt * s_next.lmd + s_next.gmm - s.gmm;
  kkt_residual.la.noalias()   += dt * s_next.gmm;
}


void linearizeUnconstrForwardEulerTerminal(const SplitSolution& s, 
                                           SplitKKTResidual& kkt_residual) {
  kkt_residual.lq().noalias() -= s.lmd;
  kkt_residual.lv().noalias() -= s.gmm;
}


void linearizeUnconstrBackwardEuler(const double dt, 
                                    const Eigen::VectorXd& q_prev, 
                                    const Eigen::VectorXd& v_prev, 
                                    const SplitSolution& s, 
                                    const SplitSolution& s_next,
                                    SplitKKTMatrix& kkt_matrix, 
                                    SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  evalUnconstrBackwardEuler(dt, q_prev, v_prev, s, kkt_residual);
  kkt_residual.lq().noalias() += s_next.lmd - s.lmd;
  kkt_residual.lv().noalias() += dt * s.lmd - s.gmm + s_next.gmm;
  kkt_residual.la.noalias()   += dt * s.gmm;
}


void linearizeUnconstrBackwardEulerTerminal(const double dt, 
                                            const Eigen::VectorXd& q_prev, 
                                            const Eigen::VectorXd& v_prev, 
                                            const SplitSolution& s, 
                                            SplitKKTMatrix& kkt_matrix, 
                                            SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  evalUnconstrBackwardEuler(dt, q_prev, v_prev, s, kkt_residual);
  kkt_residual.lq().noalias() -= s.lmd;
  kkt_residual.lv().noalias() += dt * s.lmd - s.gmm;
  kkt_residual.la.noalias()   += dt * s.gmm;
}


void evalUnconstrForwardEuler(const double dt, const SplitSolution& s, 
                              const SplitSolution& s_next, 
                              SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  kkt_residual.Fq() = s.q + dt * s.v - s_next.q;
  kkt_residual.Fv() = s.v + dt * s.a - s_next.v;
}


void evalUnconstrBackwardEuler(const double dt, const Eigen::VectorXd& q_prev, 
                               const Eigen::VectorXd& v_prev, 
                               const SplitSolution& s, 
                               SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  kkt_residual.Fq() = q_prev - s.q + dt * s.v;
  kkt_residual.Fv() = v_prev - s.v + dt * s.a;
}

} // namespace robotoc 