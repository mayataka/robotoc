#include <iostream>
#include <string>
#include <chrono>

#include "Eigen/Core"

#include "ocp/ocp.hpp"
#include "robot/robot.hpp"
#include "cost_function.hpp"
#include "constraints.hpp"


int main() {
  srand((unsigned int) time(0));
  const std::string urdf_file_name = "../crane_x7_description/urdf/crane_x7.urdf";
  idocp::Robot robot(urdf_file_name);
  const Eigen::VectorXd q_ref = Eigen::VectorXd::Constant(robot.dimq(), 0);
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Zero(robot.dimv());
  idocp::cranex7::CostFunction cost(robot, q_ref);
  idocp::cranex7::Constraints constraints(robot);
  const double T = 1;
  const unsigned int N = 50;
  const unsigned int num_proc = 1;
  idocp::OCP ocp_(robot, &cost, &constraints, T, N, num_proc);
  const double t = 0;
  Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  if (q[3] > 0.04) {
    q[3] = 0.04;
  }

  robot.printRobotModel();
  ocp_.setStateTrajectory(q, v);
  std::cout << ocp_.KKTError(t, q, v) << std::endl;
  const int num_iteration = 10;
  std::chrono::system_clock::time_point start_clock, end_clock;
  start_clock = std::chrono::system_clock::now();
  for (int i=0; i<num_iteration; ++i) {
    ocp_.solveSQP(t, q, v, true);
    std::cout << ocp_.KKTError(t, q, v) << std::endl;
  }
  end_clock = std::chrono::system_clock::now();
  // ocp_.printSolution();
  std::cout << "q = " << q.transpose() << std::endl;
  std::cout << "v = " << v.transpose() << std::endl;
  std::cout << "total CPU time: " << 1e-03 * std::chrono::duration_cast<std::chrono::microseconds>(end_clock-start_clock).count() << "[ms]" << std::endl;
  std::cout << "CPU time per update: " << 1e-03 * std::chrono::duration_cast<std::chrono::microseconds>(end_clock-start_clock).count() / num_iteration << "[ms]" << std::endl;

  return 0;
}