#include "idocp/utils/simulation_data_saver.hpp"

#include <stdexcept>

namespace idocp {

SimulationDataSaver::SimulationDataSaver(const std::string& save_dir_path, 
                                         const std::string& save_file_name) 
  : q_file_(save_dir_path+"/"+save_file_name+"_q.dat"),
    v_file_(save_dir_path+"/"+save_file_name+"_v.dat"),
    tau_file_(save_dir_path+"/"+save_file_name+"_tau.dat"),
    KKT_error_file_(save_dir_path+"/"+save_file_name+"_KKT_error.dat"),
    conditions_file_(save_dir_path+"/"+save_file_name+"conditions.dat") {
}


SimulationDataSaver::~SimulationDataSaver() {
  q_file_.close();
  v_file_.close();
  tau_file_.close();
  KKT_error_file_.close();
  conditions_file_.close();
}


void SimulationDataSaver::save(const Eigen::VectorXd& q, 
                               const Eigen::VectorXd& v, 
                               const Eigen::VectorXd& tau, 
                               const double KKT_error) {
  try {
    if (KKT_error < 0) {
      throw std::out_of_range(
          "Invalid argument: KKT error must be non negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  for (int i=0; i<q.size(); ++i) {
    q_file_ << q[i] << " ";
  }
  q_file_ << "\n";
  for (int i=0; i<v.size(); ++i) {
    v_file_ << v[i] << " ";
  }
  v_file_ << "\n";
  for (int i=0; i<tau.size(); ++i) {
    tau_file_ << tau[i] << " ";
  }
  tau_file_ << "\n";
  KKT_error_file_ << KKT_error << "\n";
}


void SimulationDataSaver::saveConditions(
    const double simulation_time_in_sec, 
    const double sampling_period_in_millisec, 
    const double CPU_time_per_update_in_millisec) {
  try {
    if (simulation_time_in_sec <= 0) {
      throw std::out_of_range(
          "Invalid argument: simulation_time_in_sec must be positive!");
    }
    if (sampling_period_in_millisec <= 0) {
      throw std::out_of_range(
          "Invalid argument: sampling_period_in_millisec must be positive!");
    }
    if (CPU_time_per_update_in_millisec <= 0) {
      throw std::out_of_range(
          "Invalid argument: CPU_time_per_update_in_millisec must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  conditions_file_ << "simulation time: " << simulation_time_in_sec << "[s]\n";
  conditions_file_ << "sampling time: " 
                  << sampling_period_in_millisec << "[ms]\n";
  conditions_file_ << "CPU time per update: " 
                  << CPU_time_per_update_in_millisec << "[ms]\n";
}

} // namespace idocp