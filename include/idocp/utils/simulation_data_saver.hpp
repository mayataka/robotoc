#ifndef IDOCP_UTILS_SIMULATION_DATA_SAVER_HPP_
#define IDOCP_UTILS_SIMULATION_DATA_SAVER_HPP_

#include <string>
#include <fstream>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"

namespace idocp {

class SimulationDataSaver {
public:
  SimulationDataSaver(const std::string& save_dir_path, 
                      const std::string& save_file_name);

  ~SimulationDataSaver();

  void save(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
            const Eigen::VectorXd& tau, const double KKT_error);

  void saveConditions(const double simulation_time_in_sec, 
                      const double sampling_period_in_millisec, 
                      const double CPU_time_per_update_in_millisec);

private:
  std::ofstream q_file_, v_file_, tau_file_, KKT_error_file_, 
                conditions_file_;

};

} // namespace idocp


#endif // IDOCP_UTILS_SIMULATION_DATA_SAVER_HPP_