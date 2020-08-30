#include <iostream>
#include <string>
#include <memory>
#include <chrono>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/parnmpc.hpp"
#include "idocp/ocp/ocp.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/contact_cost.hpp"
#include "idocp/cost/task_space_3d_cost.hpp"
#include "idocp/constraints/constraints.hpp"

#include "../utils/joint_constraints_factory.hpp"
#include "../utils/ocp_benchmarker.hpp"

namespace ocpbenchmark {
namespace cranex7 {

void BenchmarkWithoutContacts() {
  srand((unsigned int) time(0));
  const std::string urdf_file_name = "../urdf/crane_x7.urdf";
  idocp::Robot robot(urdf_file_name);
  auto cost = std::make_shared<idocp::CostFunction>();
  auto joint_cost = std::make_shared<idocp::JointSpaceCost>(robot);
  joint_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  joint_cost->set_vf_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  joint_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  joint_cost->set_u_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.00));
  const int end_effector_frame_id = 26;
  auto task_3d_cost = std::make_shared<idocp::TaskSpace3DCost>(robot, end_effector_frame_id);
  task_3d_cost->set_q_weight(Eigen::Vector3d::Constant(1000));
  task_3d_cost->set_qf_weight(Eigen::Vector3d::Constant(1000));
  Eigen::Vector3d q_ref;
  q_ref << 0, 0, 0.3;
  task_3d_cost->set_q_ref(q_ref);
  cost->push_back(joint_cost);
  cost->push_back(task_3d_cost);
  idocp::JointConstraintsFactory constraints_factory(robot);
  auto constraints = constraints_factory.create();
  const double T = 1;
  const int N = 20;
  const int num_proc = 4;
  const double t = 0;
  const Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(robot.dimv());
  idocp::OCPBenchmarker<idocp::OCP> ocp_benchmarker("OCP for task-space control of cranex7",
                                                    robot, cost, constraints, T, N, num_proc);
  ocp_benchmarker.setInitialGuessSolution(t, q, v);
  ocp_benchmarker.testConvergence(t, q, v, 30, false);
  ocp_benchmarker.printSolution();
  ocp_benchmarker.testCPUTime(t, q, v, 1000);
  idocp::OCPBenchmarker<idocp::ParNMPC> parnmpc_benchmarker("ParNMPC for task-space control of cranex7",
                                                            robot, cost, constraints, T, N, num_proc);
  parnmpc_benchmarker.setInitialGuessSolution(t, q, v);
  parnmpc_benchmarker.testConvergence(t, q, v, 30, false);
  parnmpc_benchmarker.testCPUTime(t, q, v, 1000);
}

} // namespace cranex7
} // namespace ocpbenchmark


int main() {
  ocpbenchmark::cranex7::BenchmarkWithoutContacts();
  return 0;
}