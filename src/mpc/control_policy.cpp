#include "robotoc/mpc/control_policy.hpp"

#include <stdexcept>


namespace robotoc {

ControlPolicy::ControlPolicy() 
  : t(0),
    tauJ(),
    qJ(),
    dqJ(),
    Kp(),
    Kd() {
}


ControlPolicy::ControlPolicy(const OCPSolver& ocp_solver, const double t) 
  : ControlPolicy() {
  set(ocp_solver, t);
}


void ControlPolicy::set(const OCPSolver& ocp_solver, const double _t) {
  const auto& time_discretization = ocp_solver.getTimeDiscretization();
  const auto& solution = ocp_solver.getSolution();
  const auto& lqr_policy = ocp_solver.getLQRPolicy();
  const int N = time_discretization.N();
  const int dimu = solution[0].u.size();
  t = _t;
  if (t < time_discretization[0].t) {
    tauJ = solution[0].u;
    qJ   = solution[0].q.tail(dimu);
    dqJ  = solution[0].v.tail(dimu);
    Kp   = lqr_policy[0].Kq().rightCols(dimu);
    Kd   = lqr_policy[0].Kv().rightCols(dimu);
    return;
  }
  for (int i=1; i<N-1; ++i) {
    if ((t < time_discretization[i].t) && (time_discretization[i-1].type != GridType::Impact)) {
      const double dt = time_discretization[i].t - time_discretization[i-1].t;
      const double alpha = (time_discretization[i].t - t) / dt;
      tauJ = alpha * solution[i-1].u + (1.0-alpha) * solution[i].u;
      qJ   = alpha * solution[i-1].q.tail(dimu) 
              + (1.0-alpha) * solution[i].q.tail(dimu);
      dqJ  = alpha * solution[i-1].v.tail(dimu) 
              + (1.0-alpha) * solution[i].v.tail(dimu);
      Kp   = alpha * lqr_policy[i-1].Kq().rightCols(dimu) 
              + (1.0-alpha) * lqr_policy[i].Kq().rightCols(dimu);
      Kd   = alpha * lqr_policy[i-1].Kv().rightCols(dimu) 
              + (1.0-alpha) * lqr_policy[i].Kv().rightCols(dimu);
      return;
    }
  }
  tauJ = solution[N-1].u;
  qJ   = solution[N-1].q.tail(dimu);
  dqJ  = solution[N-1].v.tail(dimu);
  Kp   = lqr_policy[N-1].Kq().rightCols(dimu);
  Kd   = lqr_policy[N-1].Kv().rightCols(dimu);
}


void ControlPolicy::disp(std::ostream& os) const {
  os << "ControlPolicy: \n";
  os << "  t:    " << t << " \n";
  os << "  tauJ: " << tauJ.transpose() << " \n";
  os << "  qJ:   " << qJ.transpose() << " \n";
  os << "  dqJ:  " << dqJ.transpose() << " \n";
  os << "  Kp: \n" << Kp << " \n";
  os << "  Kd: \n" << Kd << " \n";
}


std::ostream& operator<<(std::ostream& os, 
                          const ControlPolicy& control_policy) {
  control_policy.disp(os);
  return os;
}

} // namespace robotoc 