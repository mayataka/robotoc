#ifndef IDOCP_SIMULATOR_SIMULATOR_HPP_
#define IDOCP_SIMULATOR_SIMULATOR_HPP_

#include <string>

#include "Eigen/Core"

#include "robot/robot.hpp"
#include "ocp/mpc.hpp"

#include "runge_kutta.hpp"
#include "simulation_data_saver.hpp"


namespace idocp {
namespace simulator {

class Simulator {
public:
  Simulator(Robot& robot, const std::string& save_dir_path, 
            const std::string& save_file_name);

  void run(MPC& mpc, const double simulation_time_in_sec, 
           const double sampling_period_in_sec, 
           const double simulation_start_time_in_sec, 
           const Eigen::VectorXd& q_initial, const Eigen::VectorXd& v_initial);

private:
  RungeKutta runge_kutta_;
  SimulationDataSaver data_saver_;

};

} // namespace simulator
} // namespace idocp


#endif // IDOCP_SIMULATOR_SIMULATOR_HPP_