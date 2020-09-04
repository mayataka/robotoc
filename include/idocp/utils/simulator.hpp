#ifndef IDOCP_UTILS_SIMULATOR_HPP_
#define IDOCP_UTILS_SIMULATOR_HPP_

#include <string>

#include "Eigen/Core"

#include "idocp/ocp/mpc.hpp"
#include "idocp/robot/robot.hpp"

#include "idocp/utils/runge_kutta.hpp"
#include "idocp/utils/simulation_data_saver.hpp"


namespace idocp {

class Simulator {
public:
Simulator(Robot& robot, const std::string& save_dir_path, 
          const std::string& save_file_name);


template<typename OCPTypeDerived>
void run(MPC<OCPTypeDerived>& mpc, const double simulation_time_in_sec, 
         const double sampling_period_in_sec, 
         const double simulation_start_time_in_sec, 
         const Eigen::VectorXd& q_initial, const Eigen::VectorXd& v_initial);

private:
  RungeKutta runge_kutta_;
  SimulationDataSaver data_saver_;

};

} // namespace idocp

#include "simulator.hxx"

#endif // IDOCP_UTILS_SIMULATOR_HPP_