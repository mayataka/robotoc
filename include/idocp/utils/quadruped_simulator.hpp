#ifndef IDOCP_UTILS_QUADRUPED_SIMULATOR_HPP_
#define IDOCP_UTILS_QUADRUPED_SIMULATOR_HPP_

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

class MPCCallbackBase {
  MPCCallbackBase() {} 
  virtual ~MPCCallbackBase() = 0;

  virtual void init(const double t, const Eigen::VectorXd& q, 
                    const Eigen::VectorXd& v, 
                    OCPSolver& ocp_solver) {}

  virtual void callback(const double t, const Eigen::VectorXd& q, 
                        const Eigen::VectorXd& v, OCPSolverType& ocp_solver) {}
};


class QuadrupedSimulator {
public:
  QuadrupedSimulator(const std::string& path_to_raisim_activation_key,
                     const std::string& path_to_urdf_for_raisim, 
                     const std::string& save_dir_path, 
                     const std::string& save_file_name);

  void run(OCPSolver& ocp_solver, std::shared_ptr<MPCCallBackBase>& callback,
           const double simulation_time_in_sec, 
           const double sampling_period_in_sec, 
           const double simulation_start_time_in_sec, 
           const Eigen::VectorXd& q_initial, const Eigen::VectorXd& v_initial,
           const bool visualization=false, const bool recording=false);

  static void SetupRaisimOgre(raisim::World& world);

  static void RaiSimOgreCallback();

private:
  SimulationDataSaver data_saver_;
  std::string path_to_raisim_activation_key_, path_to_urdf_for_raisim_, 
              save_dir_path_, save_file_name_;

};

} // namespace idocp

#include "idocp/utils/quadruped_simulator.hxx"

#endif // IDOCP_UTILS_QUADRUPED_SIMULATOR_HPP_ 