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

  void runSimulation(MPC& mpc, const double t_sim, const double sampling_period,   
                     const double t_ini, const Eigen::VectorXd& q_ini, 
                     const Eigen::VectorXd& v_ini);

private:
  RungeKutta runge_kutta_;
  SimulationDataSaver data_saver_;

};

} // namespace simulator
} // namespace idocp


#endif // IDOCP_SIMULATOR_SIMULATOR_HPP_