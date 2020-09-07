#ifndef IDOCP_QUADRUPED_UTILS_SIMULATOR_HPP_
#define IDOCP_QUADRUPED_UTILS_SIMULATOR_HPP_

#include <string>

#include "Eigen/Core"

#include "idocp/ocp/mpc.hpp"
#include "idocp/robot/robot.hpp"

#include "idocp/utils/runge_kutta.hpp"
#include "idocp/utils/simulation_data_saver.hpp"

#include "raisim/World.hpp"
#include "raisim/OgreVis.hpp"
#include "raisim/helper.hpp"


namespace idocp {

class QuadrupedSimulator {
public:
QuadrupedSimulator(const std::string& path_to_urdf_for_raisim, 
                   const std::string& save_dir_path, 
                   const std::string& save_file_name);

template<typename OCPTypeDerived>
void run(MPC<OCPTypeDerived>& mpc, const double simulation_time_in_sec, 
         const double sampling_period_in_sec, 
         const double simulation_start_time_in_sec, 
         const Eigen::VectorXd& q_initial, const Eigen::VectorXd& v_initial);

template<typename OCPTypeDerived>
void viz(MPC<OCPTypeDerived>& mpc, const double simulation_time_in_sec, 
         const double sampling_period_in_sec, 
         const double simulation_start_time_in_sec, 
         const Eigen::VectorXd& q_initial, const Eigen::VectorXd& v_initial,
         const bool recording=false);

static void SetupCallback();

private:
  SimulationDataSaver data_saver_;
  raisim::World raisim_world_;
  raisim::ArticulatedSystem* raisim_robot_;
  raisim::Ground* raisim_ground_;
  std::string save_dir_path_, save_file_name_;
};

} // namespace idocp

#include "quadruped_simulator.hxx"

#endif // IDOCP_QUADRUPED_UTILS_SIMULATOR_HPP_