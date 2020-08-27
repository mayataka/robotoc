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
  joint_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  joint_cost->set_q_ref(Eigen::VectorXd::Constant(robot.dimv(), 2));
  joint_cost->set_qf_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  joint_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  joint_cost->set_vf_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  joint_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  joint_cost->set_u_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  cost->push_back(joint_cost);
  idocp::JointConstraintsFactory constraints_factory(robot);
  auto constraints = constraints_factory.create();
  const double T = 1;
  const int N = 20;
  const int num_proc = 4;
  const double t = 0;
  const Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(robot.dimv());
  idocp::OCPBenchmarker<idocp::OCP> ocp_benchmarker("OCP for cranex7 without contacts",
                                                    robot, cost, constraints, T, N, num_proc);
  ocp_benchmarker.setInitialGuessSolution(t, q, v);
  ocp_benchmarker.testConvergence(t, q, v, 20, false);
  ocp_benchmarker.testCPUTime(t, q, v, 1000);
  idocp::OCPBenchmarker<idocp::ParNMPC> parnmpc_benchmarker("ParNMPC for cranex7 without contacts",
                                                            robot, cost, constraints, T, N, num_proc);
  parnmpc_benchmarker.setInitialGuessSolution(t, q, v);
  parnmpc_benchmarker.testConvergence(t, q, v, 20, false);
  parnmpc_benchmarker.testCPUTime(t, q, v, 1000);
}

} // namespace cranex7
} // namespace ocpbenchmark


int main() {
  ocpbenchmark::cranex7::BenchmarkWithoutContacts();
  return 0;
}