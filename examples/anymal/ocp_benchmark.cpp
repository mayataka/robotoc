#include <iostream>
#include <string>
#include <memory>
#include <chrono>

#include "Eigen/Core"

#include "idocp/ocp/ocp.hpp"
#include "idocp/robot/robot.hpp"
#include "idocp/manipulator/cost_function.hpp"
#include "idocp/manipulator/constraints.hpp"


namespace ocpbenchmark {
namespace anymal {

void CPUTime_with_contacts() {
  srand((unsigned int) time(0));
  std::vector<int> contact_frames = {14, 24, 34, 44};
  const double baumgarte_weight_on_velocity = 10;
  const double baumgarte_weight_on_position = 10;
  const std::string urdf_file_name = "../anymal/anymal.urdf";
  idocp::Robot robot(urdf_file_name, contact_frames, 
                     baumgarte_weight_on_velocity, 
                     baumgarte_weight_on_position);
  std::shared_ptr<idocp::CostFunctionInterface> cost 
      = std::make_shared<idocp::manipulator::CostFunction>(robot);
  std::shared_ptr<idocp::ConstraintsInterface> constraints
      = std::make_shared<idocp::manipulator::Constraints>(robot);
  const double T = 1;
  const unsigned int N = 50;
  const unsigned int num_proc = 4;
  idocp::OCP ocp_(robot, cost, constraints, T, N, num_proc);
  const double t = 0;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  q << 0, 0, 4.8, 0, 0, 0, 1, 
       0.0315, 0.4, -0.8, 
       0.0315, -0.4, 0.8, 
       -0.0315, 0.4, -0.8,
       -0.0315, -0.4, 0.8;
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(robot.dimv());
  ocp_.setStateTrajectory(q, v);
  std::vector<bool> contact_status = {true, true, true, true};
  std::vector<std::vector<bool>> contact_sequence = {N, contact_status};
  ocp_.setContactSequence(contact_sequence);
  const int num_iteration = 50;
  std::chrono::system_clock::time_point start_clock, end_clock;
  start_clock = std::chrono::system_clock::now();
  for (int i=0; i<num_iteration; ++i) {
    ocp_.solveLQR(t, q, v, false);
  }
  end_clock = std::chrono::system_clock::now();
  std::cout << "---------- OCP benchmark ----------" << std::endl;
  std::cout << "model: anymal" << std::endl;
  std::cout << "dimq = " << robot.dimq() << std::endl;
  std::cout << "dimv = " << robot.dimv() << std::endl;
  std::cout << "max_dimf = " << robot.max_dimf() << std::endl;
  std::cout << "N = " << N << std::endl;
  std::cout << "T = " << T << std::endl;
  std::cout << "number of threads = " << num_proc << std::endl;
  std::cout << "total CPU time: " << 1e-03 * std::chrono::duration_cast<std::chrono::microseconds>(end_clock-start_clock).count() << "[ms]" << std::endl;
  std::cout << "CPU time per update: " << 1e-03 * std::chrono::duration_cast<std::chrono::microseconds>(end_clock-start_clock).count() / num_iteration << "[ms]" << std::endl;
  std::cout << "-----------------------------------" << std::endl;
  std::cout << std::endl;
}


void KKTError_with_contacts() {
  srand((unsigned int) time(0));
  std::vector<int> contact_frames = {14, 24, 34, 44};
  const double baumgarte_weight_on_velocity = 10;
  const double baumgarte_weight_on_position = 10;
  const std::string urdf_file_name = "../anymal/anymal.urdf";
  idocp::Robot robot(urdf_file_name, contact_frames, 
                     baumgarte_weight_on_velocity, 
                     baumgarte_weight_on_position);
  std::shared_ptr<idocp::CostFunctionInterface> cost 
      = std::make_shared<idocp::manipulator::CostFunction>(robot);
  std::shared_ptr<idocp::ConstraintsInterface> constraints
      = std::make_shared<idocp::manipulator::Constraints>(robot);
  const double T = 1;
  const unsigned int N = 50;
  const unsigned int num_proc = 4;
  const double t = 0;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  q << 0, 0, 4.8, 0, 0, 0, 1, 
       0.0315, 0.4, -0.8, 
       0.0315, -0.4, 0.8, 
       -0.0315, 0.4, -0.8,
       -0.0315, -0.4, 0.8;
  robot.normalizeConfiguration(q);
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(robot.dimv());
  robot.updateKinematics(q, v, Eigen::VectorXd::Zero(robot.dimv()));
  robot.setContactPointsByCurrentKinematics();
  idocp::OCP ocp_(robot, cost, constraints, T, N, num_proc);
  std::vector<bool> contact_status = {true, true, true, true};
  std::vector<std::vector<bool>> contact_sequence = {N, contact_status};
  ocp_.setContactSequence(contact_sequence);
  ocp_.setStateTrajectory(q, v);
  std::cout << "---------- OCP benchmark ----------" << std::endl;
  std::cout << "model: anymal" << std::endl;
  std::cout << "dimq = " << robot.dimq() << std::endl;
  std::cout << "dimv = " << robot.dimv() << std::endl;
  std::cout << "max_dimf = " << robot.max_dimf() << std::endl;
  std::cout << "q = " << q.transpose() << std::endl;
  std::cout << "v = " << v.transpose() << std::endl;
  std::cout << "Initial KKT error = " << ocp_.KKTError(t, q, v) << std::endl;
  const int num_iteration = 50;
  for (int i=0; i<num_iteration; ++i) {
    ocp_.solveLQR(t, q, v, true);
    std::cout << "KKT error at iteration " << i << " = " << ocp_.KKTError(t, q, v) << std::endl;
  }
  std::cout << "-----------------------------------" << std::endl;
  std::cout << std::endl;
}

} // namespace anymal
} // namespace ocpbenchmark


int main() {
  ocpbenchmark::anymal::CPUTime_with_contacts();
  ocpbenchmark::anymal::KKTError_with_contacts();
  return 0;
}
