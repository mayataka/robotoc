#include <iostream>
#include <string>

#include "Eigen/Core"

#include "ocp/mpc.hpp"
#include "robot/robot.hpp"
#include "manipulator/cost_function.hpp"
#include "manipulator/constraints.hpp"

#include "mpc_simulator/simulator.hpp"


int main() {
  srand((unsigned int) time(0));
  const std::string urdf_file_name = "../urdf/iiwa14.urdf";
  idocp::Robot robot(urdf_file_name);
  idocp::manipulator::CostFunction cost(robot);
  idocp::manipulator::Constraints constraints(robot);
  Eigen::VectorXd q_ref = 1.5 * Eigen::VectorXd::Random(robot.dimq());
  cost.set_q_ref(q_ref);
  const double T = 1;
  const unsigned int N = 18;
  const unsigned int num_proc = 4;
  idocp::MPC mpc(robot, &cost, &constraints, T, N, num_proc);
  const double t = 0;
  Eigen::VectorXd q0 = Eigen::VectorXd::Zero(robot.dimq());
  Eigen::VectorXd v0 = Eigen::VectorXd::Zero(robot.dimv());
  mpc.initializeSolution(t, q0, v0, 0);

  const double t0 = 0;
  const double simulation_time = 10;
  const double sampling_period = 0.001;
  const std::string save_dir_path = "../sim_result";
  const std::string save_file_name = "iiwa14";
  idocp::simulator::Simulator simulator(robot, save_dir_path, save_file_name);
  simulator.run(mpc, simulation_time, sampling_period, t0, q0, v0);
  return 0;
}