#include <iostream>
#include <string>
#include <memory>
#include <chrono>

#include "Eigen/Core"

#include "idocp/ocp/ocp.hpp"
#include "idocp/robot/robot.hpp"
#include "idocp/manipulator/cost_function.hpp"
#include "idocp/manipulator/constraints.hpp"


int main() {
  srand((unsigned int) time(0));
  std::vector<int> contact_frames(1);
  contact_frames = {18};
  const double baumgarte_weight_on_velocity = 10;
  const double baumgarte_weight_on_position = 100;
  const std::string urdf_file_name = "../urdf/iiwa14.urdf";
  idocp::Robot robot(urdf_file_name, contact_frames, 
                     baumgarte_weight_on_velocity, 
                     baumgarte_weight_on_position);
  std::shared_ptr<idocp::CostFunctionInterface> cost 
      = std::make_shared<idocp::manipulator::CostFunction>(robot);
  std::shared_ptr<idocp::ConstraintsInterface> constraints
      = std::make_shared<idocp::manipulator::Constraints>(robot);
  Eigen::VectorXd q_ref = 1.5 * Eigen::VectorXd::Random(robot.dimq());
  const double T = 1;
  const unsigned int N = 50;
  const unsigned int num_proc = 1;
  idocp::OCP ocp_(robot, cost, constraints, T, N, num_proc);
  const double t = 0;
  Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  ocp_.setStateTrajectory(q, v);
  std::cout << ocp_.KKTError(t, q, v) << std::endl;
  const int num_iteration = 10;
  std::chrono::system_clock::time_point start_clock, end_clock;
  start_clock = std::chrono::system_clock::now();
  for (int i=0; i<num_iteration; ++i) {
    ocp_.solveLQR(t, q, v, false);
    std::cout << ocp_.KKTError(t, q, v) << std::endl;
  }
  end_clock = std::chrono::system_clock::now();
  // ocp_.printSolution();
  std::cout << "q = " << q.transpose() << std::endl;
  std::cout << "v = " << v.transpose() << std::endl;
  std::cout << "q_ref = " << q_ref.transpose() << std::endl;
  std::cout << "total CPU time: " << 1e-03 * std::chrono::duration_cast<std::chrono::microseconds>(end_clock-start_clock).count() << "[ms]" << std::endl;
  std::cout << "CPU time per update: " << 1e-03 * std::chrono::duration_cast<std::chrono::microseconds>(end_clock-start_clock).count() / num_iteration << "[ms]" << std::endl;

  return 0;
}
