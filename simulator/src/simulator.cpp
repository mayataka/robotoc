#include "simulator.hpp"


namespace idocp {
namespace simulator {

Simulator::Simulator(Robot& robot, const std::string& save_dir_path, 
                     const std::string& save_file_name)
  : runge_kutta_(robot),
    data_saver_(save_dir_path, save_file_name) {
}


void Simulator::runSimulation(MPC& mpc, const double t_sim, 
                              const double sampling_period, const double t_ini, 
                              const Eigen::VectorXd& q_ini, 
                              const Eigen::VectorXd& v_ini) {
  double t = t_ini;
  Eigen::VectorXd q = q_ini;
  Eigen::VectorXd q_next = q;
  Eigen::VectorXd v = v_ini;
  Eigen::VectorXd v_next = v;
  Eigen::VectorXd u = Eigen::VectorXd::Zero(v_ini.size());
  mpc.getControlInput(u);
  while (t < t_sim) {
    data_saver_.save(q, v, u, mpc.KKTError(t, q, v));
    runge_kutta_.integrate(sampling_period, q, v, u, q_next, v_next);
    mpc.updateSolution(t, q, v);
    mpc.getControlInput(u);
    t += sampling_period;
    q = q_next;
    v = v_next;
  }
}

} // namespace simulator
} // namespace idocp