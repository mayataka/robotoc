#ifndef IDOCP_UTILS_MPC_SIMULATOR_HXX_
#define IDOCP_UTILS_MPC_SIMULATOR_HXX_

#include "idocp/utils/mpc_simulator.hpp"

#include <chrono>
#include <stdexcept>

namespace idocp {

inline MPCSimulator::MPCSimulator(Robot& robot, 
                                  const std::string& path_to_save_dir, 
                                  const std::string& save_file_name)
  : runge_kutta_(robot),
    data_saver_(path_to_save_dir, save_file_name) {
}


template<typename OCPTypeDerived>
inline void MPCSimulator::run(MPC<OCPTypeDerived>& mpc, 
                              const double simulation_time_in_sec, 
                              const double sampling_period_in_sec, 
                              const double simulation_start_time_in_sec, 
                              const Eigen::VectorXd& q_initial, 
                              const Eigen::VectorXd& v_initial) {
  try {
    if (simulation_time_in_sec <= 0) {
      throw std::out_of_range(
          "Invalid argument: simulation_time_in_sec must be positive!");
    }
    if (sampling_period_in_sec <= 0) {
      throw std::out_of_range(
          "Invalid argument: sampling_period_in_sec must be positive!");
    }
    if (simulation_start_time_in_sec < 0) {
      throw std::out_of_range(
          "Invalid argument: simulation_start_time_in_sec must be non negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  Eigen::VectorXd q = q_initial;
  Eigen::VectorXd q_next = q;
  Eigen::VectorXd v = v_initial;
  Eigen::VectorXd v_next = v;
  Eigen::VectorXd u = Eigen::VectorXd::Zero(v_initial.size());
  mpc.getControlInput(u);
  std::chrono::system_clock::time_point start_clock, end_clock;
  double CPU_time_total_in_sec = 0;
  for (double t=0; t<simulation_time_in_sec; t+=sampling_period_in_sec) {
		mpc.computeKKTResidual(t, q, v);
    data_saver_.save(q, v, u, mpc.KKTError());
    runge_kutta_.integrate(sampling_period_in_sec, q, v, u, q_next, v_next);
    start_clock = std::chrono::system_clock::now();
    mpc.updateSolution(t, q, v);
    mpc.getControlInput(u);
    end_clock = std::chrono::system_clock::now();
    const double CPU_time_in_sec
        = 1.0e-06 * std::chrono::duration_cast<std::chrono::microseconds>(
            end_clock-start_clock).count();
    CPU_time_total_in_sec += CPU_time_in_sec;
    q = q_next;
    v = v_next;
  }
  const int num_simulation_update 
      = (int)(simulation_time_in_sec/sampling_period_in_sec);
  const double CPU_time_per_update_in_sec 
      = CPU_time_total_in_sec / num_simulation_update;
  data_saver_.saveConditions(simulation_time_in_sec, 
                             1000*sampling_period_in_sec,
                             1000*CPU_time_per_update_in_sec);
  std::cout << "simulation time: " << simulation_time_in_sec << "[s]" 
            << std::endl;
  std::cout << "sampling time: " << 1000 * sampling_period_in_sec << "[ms]" 
            << std::endl;
  std::cout << "CPU time per update: " << 1000 * CPU_time_per_update_in_sec 
            << "[ms]" << std::endl;
}

} // namespace idocp

#endif // IDOCP_UTILS_MPC_SIMULATOR_HXX_ 