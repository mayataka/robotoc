#include <iostream>
#include <string>

#include "Eigen/Core"

#include "ocp/mpc.hpp"
#include "robot/robot.hpp"
#include "cost_function.hpp"
#include "constraints.hpp"

#include "simulator.hpp"


int main() {
  srand((unsigned int) time(0));
  const std::string urdf_file_name = "../urdf/iiwa14.urdf";
  idocp::Robot robot(urdf_file_name);
  const Eigen::VectorXd q_ref = Eigen::VectorXd::Constant(robot.dimq(), -2);
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Zero(robot.dimv());
  idocp::iiwa14::CostFunction cost(robot, q_ref);
  idocp::iiwa14::Constraints constraints(robot);
  const double T = 1;
  const unsigned int N = 50;
  const unsigned int num_proc = 1;
  idocp::MPC mpc(robot, &cost, &constraints, T, N, num_proc);
  const double t = 0;
  Eigen::VectorXd q0 = Eigen::VectorXd::Random(robot.dimq());
  Eigen::VectorXd v0 = Eigen::VectorXd::Random(robot.dimv());
  q0.fill(2);
  v0.fill(0);
  // mpc.initializeSolution(t, q1, v0);

  const double t0 = 0;
  const double t_sim = 10;
  const double sampling_period = 0.001;
  const std::string save_dir_path = "../sim_result";
  const std::string save_file_name = "iiwa14";
  idocp::simulator::Simulator simulator(robot, save_dir_path, save_file_name);
  simulator.runSimulation(mpc, t_sim, sampling_period, t0, q0, v0);
  return 0;
}