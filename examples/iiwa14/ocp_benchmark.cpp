#include <iostream>
#include <string>
#include <memory>
#include <chrono>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
// #include "idocp/ocp/parnmpc.hpp"
#include "idocp/ocp/ocp.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/contact_force_cost.hpp"
#include "idocp/constraints/constraints.hpp"

#include "idocp/utils/joint_constraints_factory.hpp"
#include "idocp/utils/ocp_benchmarker.hpp"

namespace ocpbenchmark {
namespace iiwa14 {

void BenchmarkWithoutContacts() {
  srand((unsigned int) time(0));
  const std::string urdf_file_name = "../urdf/iiwa14.urdf";
  idocp::Robot robot(urdf_file_name);
  auto cost = std::make_shared<idocp::CostFunction>();
  auto joint_cost = std::make_shared<idocp::JointSpaceCost>(robot);
  joint_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  joint_cost->set_qf_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  joint_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  joint_cost->set_vf_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  joint_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  joint_cost->set_u_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.0));
  cost->push_back(joint_cost);
  idocp::JointConstraintsFactory constraints_factory(robot);
  auto constraints = constraints_factory.create();
  const double T = 1;
  const int N = 50;
  const int num_proc = 4;
  const double t = 0;
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  idocp::OCPBenchmarker<idocp::OCP> ocp_benchmarker("OCP for iiwa14 without contacts",
                                                    robot, cost, constraints, T, N, num_proc);
  ocp_benchmarker.setInitialGuessSolution(t, q, v);
  ocp_benchmarker.testConvergence(t, q, v, 20, false);
  // ocp_benchmarker.testCPUTime(t, q, v, 10000);
  // idocp::OCPBenchmarker<idocp::ParNMPC> parnmpc_benchmarker("ParNMPC for iiwa14 without contacts",
  //                                                           robot, cost, constraints, T, N, num_proc);
  // parnmpc_benchmarker.setInitialGuessSolution(t, q, v);
  // parnmpc_benchmarker.testConvergence(t, q, v, 20, false);
  // parnmpc_benchmarker.testCPUTime(t, q, v, 10000);
}


void BenchmarkWithContacts() {
  srand((unsigned int) time(0));
  std::vector<int> contact_frames = {18};
  const std::string urdf_file_name = "../urdf/iiwa14.urdf";
  idocp::Robot robot(urdf_file_name, contact_frames);
  auto cost = std::make_shared<idocp::CostFunction>();
  auto joint_cost = std::make_shared<idocp::JointSpaceCost>(robot);
  joint_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  joint_cost->set_qf_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  joint_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  joint_cost->set_vf_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  joint_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  joint_cost->set_u_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.0));
  auto contact_force_cost = std::make_shared<idocp::ContactForceCost>(robot);
  std::vector<Eigen::Vector3d> f_weight;
  for (int i=0; i<contact_frames.size(); ++i) {
    f_weight.push_back(Eigen::Vector3d::Constant(0.001));
  }
  contact_force_cost->set_f_weight(f_weight);
  cost->push_back(joint_cost);
  cost->push_back(contact_force_cost);
  idocp::JointConstraintsFactory constraints_factory(robot);
  auto constraints = constraints_factory.create();
  const double T = 1;
  const int N = 50;
  const int num_proc = 4;
  const double t = 0;
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(robot.dimv());
  joint_cost->set_q_ref(q);
  robot.updateKinematics(q, v);
  robot.setContactPointsByCurrentKinematics();
  idocp::OCPBenchmarker<idocp::OCP> ocp_benchmarker("OCP for iiwa14 with contacts",
                                                    robot, cost, constraints, T, N, num_proc);
  ocp_benchmarker.setInitialGuessSolution(t, q, v);
  ocp_benchmarker.getSolverHandle()->activateContacts({0}, 0, N);
  ocp_benchmarker.testConvergence(t, q, v, 20, false);

  const int N2 = 100;
  idocp::OCPBenchmarker<idocp::OCP> ocp_benchmarker2("OCP for iiwa14 with contacts2",
                                                    robot, cost, constraints, T, N2, num_proc);
  ocp_benchmarker2.setInitialGuessSolution(t, q, v);
  ocp_benchmarker2.getSolverHandle()->activateContacts({0}, 0, N2);
  ocp_benchmarker2.testConvergence(t, q, v, 20, false);
  // ocp_benchmarker.testCPUTime(t, q, v);
  // idocp::OCPBenchmarker<idocp::ParNMPC> parnmpc_benchmarker("ParNMPC for iiwa14 with contacts",
  //                                                           robot, cost, constraints, T, N, num_proc);
  // parnmpc_benchmarker.setInitialGuessSolution(t, q, v);
  // parnmpc_benchmarker.getSolverHandle()->activateContacts({0}, 0, N);
  // parnmpc_benchmarker.testConvergence(t, q, v, 30, true);
  // parnmpc_benchmarker.testCPUTime(t, q, v);
}


} // namespace iiwa14
} // namespace ocpbenchmark


int main() {
  ocpbenchmark::iiwa14::BenchmarkWithoutContacts();
  ocpbenchmark::iiwa14::BenchmarkWithContacts();
  return 0;
}
