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
  const std::string urdf_file_name = "../crane_x7_description/urdf/crane_x7.urdf";
  idocp::Robot robot(urdf_file_name);
  const Eigen::VectorXd q_ref = Eigen::VectorXd::Constant(robot.dimq(), 0);
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Zero(robot.dimv());
  idocp::cranex7::CostFunction cost(robot, q_ref);
  idocp::cranex7::Constraints constraints(robot);
  const double T = 1;
  const unsigned int N = 18;
  const unsigned int num_proc = 4;
  idocp::MPC mpc(robot, &cost, &constraints, T, N, num_proc);
  const double t = 0;
  Eigen::VectorXd q0 = Eigen::VectorXd::Random(robot.dimq());
  Eigen::VectorXd v0 = Eigen::VectorXd::Random(robot.dimv());
  if (q0[3] > 0.04) {
    q0[3] = 0.04;
  }
  mpc.initializeSolution(t, q0, v0, 0);

  const double t0 = 0;
  const double simulation_time = 10;
  const double sampling_period = 0.001;
  const std::string save_dir_path = "../sim_result";
  const std::string save_file_name = "cranex7";
  idocp::simulator::Simulator simulator(robot, save_dir_path, save_file_name);
  simulator.run(mpc, simulation_time, sampling_period, t0, q0, v0);
  return 0;
}